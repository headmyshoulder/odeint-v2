#ifndef PID_STEP_ADJUSTER_HPP_INCLUDED
#define PID_STEP_ADJUSTER_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>
#include <boost/numeric/odeint/stepper/detail/pid_step_adjuster_coefficients.hpp>

#include <math.h>

namespace boost{
namespace numeric{
namespace odeint{
namespace detail{

template<
size_t Steps,
class State,
class Value = double,
class Time = double,
size_t Type = H211PI
>
struct pid_step_adjuster
{
public:
    const static size_t steps = Steps;
    static double threshold() { return 0.9; };

    typedef State state_type;
    typedef Value value_type;
    typedef Time time_type;

    typedef rotating_buffer<state_type, 3> error_storage_type;
    typedef rotating_buffer<time_type, 3> time_storage_type;
    typedef pid_step_adjuster_coefficients<Type> coeff_type;

    pid_step_adjuster(double tol = 1e-6, time_type dtmax = 1.0)
    :m_dtmax(dtmax), m_error_storage(), m_time_storage(), m_init(0), m_tol(tol)
    {};

    time_type adjust_stepsize(time_type dt, const state_type &err)
    {
        m_error_storage[0] = err;
        m_time_storage[0] = dt;

        value_type ratio = 100;
        value_type r;

        for(size_t i=0; i<m_error_storage[0].size(); ++i)
        {
            if(m_init >= 2)
            {
                r = adapted_pow(fabs(m_tol/m_error_storage[0][i]), m_coeff[0]/(steps + 1)) *
                    adapted_pow(fabs(m_tol/m_error_storage[1][i]), m_coeff[1]/(steps + 1)) *
                    adapted_pow(fabs(m_tol/m_error_storage[2][i]), m_coeff[2]/(steps + 1)) *
                    adapted_pow(fabs(m_time_storage[0]/m_time_storage[1]), -m_coeff[3]/(steps + 1))*
                    adapted_pow(fabs(m_time_storage[1]/m_time_storage[2]), -m_coeff[4]/(steps + 1));
            }
            else
            {
                r = pow(fabs(m_tol/m_error_storage[0][i]), 0.7/(steps + 1)); // purely integrating controller for startup
            }

            if(r<ratio)
                ratio = r;
        }

        value_type kappa = 1.0;
        ratio = 1.0 + kappa*atan((ratio - 1)/kappa);

        if(ratio*dt >= m_dtmax)
        {
            dt = m_dtmax;
        }

        if(ratio >= threshold() )
        {
            m_error_storage.rotate();
            m_time_storage.rotate();

            ++m_init;
        }
        else
        {
            m_init = 0;
        }

        return dt * static_cast<time_type> (ratio);
    };

private:
    template<class T>
    inline value_type adapted_pow(T base, double exp)
    {
        if (exp > 0)
            return pow(base, exp);
        else
            return 1/pow(base, -exp);
    };

    time_type m_dtmax;
    error_storage_type m_error_storage;
    time_storage_type m_time_storage;

    size_t m_init;
    double m_tol;

    coeff_type m_coeff;
};

}
}
}
}
#endif