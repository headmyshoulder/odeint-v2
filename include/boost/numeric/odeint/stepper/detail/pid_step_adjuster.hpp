#ifndef PID_STEP_ADJUSTER_HPP_INCLUDED
#define PID_STEP_ADJUSTER_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>

#include <math.h>

namespace boost{
namespace numeric{
namespace odeint{
namespace detail{
template<size_t order, class State, class Time, size_t tol = 10000>
struct pid_step_adjuster
{
	public:
		typedef State state_type;
		typedef Time time_type;
		typedef rotating_buffer<state_type, 3> error_storage_type;

		pid_step_adjuster(time_type dtmax = 0.5)
		:m_dtmax(dtmax), m_error_storage(), init(0)
		{};

		double adjust_stepsize(const state_type &err, const time_type &dt)
		{
			// basic controller
			m_error_storage[0] = err;

			if(init < order-1)
			{
				++init;
				m_error_storage.rotate();
				return 1.0;
			}

			double theta = 1.0/tol;
			double ratio = 100;
			double r;

			for(size_t i=0; i<m_error_storage[0].size(); ++i)
			{
				r = pow(fabs(theta/m_error_storage[0][i]), 1.0/(order+1));
				if(r<ratio)
					ratio = r;
			}

			double kappa = 1.0; // against too big scaling steps
			ratio = 1 + kappa*atan((ratio - 1)/kappa);

			if(ratio*dt >= m_dtmax)
			{
				ratio = m_dtmax / dt;
			}

			if(ratio > 0.9)
			{
				m_error_storage.rotate();
			}
			
			return ratio;
		};

	private:
		time_type m_dtmax;
		error_storage_type m_error_storage;
		size_t init;
};

}
}
}
}
#endif