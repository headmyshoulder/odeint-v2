#ifndef PID_STEP_ADJUSTER_HPP_INCLUDED
#define PID_STEP_ADJUSTER_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>
#include <boost/numeric/odeint/stepper/detail/pid_step_adjuster_coefficients.hpp>

#include <math.h>

namespace boost{
namespace numeric{
namespace odeint{
namespace detail{

template<size_t Steps, class State, class Time, size_t Type = H211PI>
struct pid_step_adjuster
{
	public:
		typedef State state_type;
		typedef Time time_type;

		typedef rotating_buffer<state_type, 3> error_storage_type;
		typedef rotating_buffer<time_type, 3> time_storage_type;
		typedef pid_step_adjuster_coefficients<Type> coeff_type;

		pid_step_adjuster(double tol = 1e-5, time_type dtmax = 1.0)
		:m_dtmax(dtmax), m_error_storage(), m_time_storage(), m_tol(tol), m_failed(0), m_init(0)
		{};

		double adjust_stepsize(const state_type &err, const time_type &dt)
		{
			// basic controller
			m_error_storage[0] = err;
			m_time_storage[0] = dt;

			double ratio = 100;
			double r;

			for(size_t i=0; i<m_error_storage[0].size(); ++i)
			{
				if(m_init >= 2)
				{
					r = pow(fabs(m_tol/m_error_storage[0][i]), m_coeff[0]/(Steps + 1)) *
						pow(fabs(m_tol/m_error_storage[1][i]), m_coeff[1]/(Steps + 1)) *
						pow(fabs(m_tol/m_error_storage[2][i]), m_coeff[2]/(Steps + 1)) *
						pow(fabs(m_time_storage[0]/m_time_storage[1]), -m_coeff[3]/(Steps + 1))*
						pow(fabs(m_time_storage[1]/m_time_storage[2]), -m_coeff[4]/(Steps + 1));
				 }
				 else
				 {
				 	r = pow(fabs(m_tol/m_error_storage[0][i]), 0.7/(Steps + 1)); // purely integrating controller for startup
			 	}

				if(r<ratio)
					ratio = r;
			}

			double kappa = 1.0;
			ratio = 1.0 + kappa*atan((ratio - 1)/kappa);

			if(ratio*dt >= m_dtmax)
			{
				ratio = m_dtmax / dt;
			}

			if(ratio >= 0.9)
			{
				m_error_storage.rotate();
				m_time_storage.rotate();

				++m_init;
			}
			else
			{
				m_init = 0;
				m_failed ++;
				//std::cout << m_failed << std::endl;
			}

			return ratio;
		};

	private:
		time_type m_dtmax;
		error_storage_type m_error_storage;
		time_storage_type m_time_storage;

		size_t m_init;
		double m_tol;

		size_t m_failed;

		coeff_type m_coeff;
};

}
}
}
}
#endif
