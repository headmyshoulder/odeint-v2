#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_PID_STEP_ADJUSTER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_PID_STEP_ADJUSTER_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>
#include <boost/numeric/odeint/stepper/detail/pid_step_adjuster_coefficients.hpp>

#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>

#include <math.h>

namespace boost{
namespace numeric{
namespace odeint{
namespace detail{

template<
size_t Steps,
class Value = double,
class Time = double
>
struct op
{
public:
    static const size_t steps = Steps;

    typedef Value value_type;
    typedef Time time_type;

    const double beta1;
    const double beta2;
    const double beta3;
    const double alpha1;
    const double alpha2;

    const time_type t1;
    const time_type t2;
    const time_type t3;

    const double m_tol;

    op(const double tol, const double _t1, const double _t2, const double _t3,
        const double b1 = 1, const double b2 = 0, const double b3 = 0, const double a1 = 0, const double a2 = 0)
    :beta1(b1), beta2(b2), beta3(b3), alpha1(a1), alpha2(a2),
    t1(_t1), t2(_t2), t3(_t3), m_tol(tol)
    {};

    template<class T1, class T2, class T3, class T4>
    void operator()(T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4)
    {
        t1 = adapted_pow(fabs(m_tol/t2), beta1/(steps + 1)) *
            adapted_pow(fabs(m_tol/t3), beta2/(steps + 1)) *
            adapted_pow(fabs(m_tol/t4), beta3/(steps + 1)) *
            adapted_pow(fabs(t1/t2), -alpha1/(steps + 1))*
            adapted_pow(fabs(t2/t3), -alpha2/(steps + 1));
    };

    template<class T1, class T2>
    void operator()(T1 &t1, const T2 &t2)
    {
        t1 = adapted_pow(fabs(m_tol/t2), beta1/(steps + 1));
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
};

template<
size_t Steps,
class State,
class Value = double,
class Time = double,
class Algebra = typename algebra_dispatcher< State >::algebra_type,
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
    typedef Algebra algebra_type;

    typedef rotating_buffer<state_type, 3> error_storage_type;
    typedef rotating_buffer<time_type, 3> time_storage_type;
    typedef pid_step_adjuster_coefficients<Type> coeff_type;

    pid_step_adjuster(double tol = 1e-6, time_type dtmax = 1.0)
    :m_dtmax(dtmax), m_error_storage(), m_dt_storage(), m_init(0), m_tol(tol)
    {};

    time_type adjust_stepsize(time_type dt, state_type &err)
    {
        m_error_storage[0] = err;
        m_dt_storage[0] = dt;

        value_type ratio;
        // value_type r;

        if(m_init >= 2)
        {
            m_algebra.for_each4(err, m_error_storage[0], m_error_storage[1], m_error_storage[2],
                op<steps>(m_tol, m_dt_storage[0], m_dt_storage[1], m_dt_storage[2],
                    m_coeff[0], m_coeff[1], m_coeff[2], m_coeff[3], m_coeff[4]));
        }
        else
        {
            m_algebra.for_each2(err, m_error_storage[0],
                op<steps>(m_tol, m_dt_storage[0], m_dt_storage[1], m_dt_storage[2], 0.7));
        }

        ratio = m_algebra.norm_inf(err);

        /*for(size_t i=0; i<m_error_storage[0].size(); ++i)
        {
            if(m_init >= 2)
            {
                r = adapted_pow(fabs(m_tol/m_error_storage[0][i]), m_coeff[0]/(steps + 1)) *
                    adapted_pow(fabs(m_tol/m_error_storage[1][i]), m_coeff[1]/(steps + 1)) *
                    adapted_pow(fabs(m_tol/m_error_storage[2][i]), m_coeff[2]/(steps + 1)) *
                    adapted_pow(fabs(m_dt_storage[0]/m_dt_storage[1]), -m_coeff[3]/(steps + 1))*
                    adapted_pow(fabs(m_dt_storage[1]/m_dt_storage[2]), -m_coeff[4]/(steps + 1));
            }
            else
            {
                r = pow(fabs(m_tol/m_error_storage[0][i]), 0.7/(steps + 1)); // purely integrating controller for startup
            }

            if(r<ratio)
                ratio = r;
        }*/

        value_type kappa = 1.0;
        ratio = 1.0 + kappa*atan((ratio - 1)/kappa);

        if(ratio*dt >= m_dtmax)
        {
            ratio = m_dtmax / dt;
        }

        if(ratio >= threshold() )
        {
            m_error_storage.rotate();
            m_dt_storage.rotate();

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

    algebra_type m_algebra;

    time_type m_dtmax;
    error_storage_type m_error_storage;
    time_storage_type m_dt_storage;

    size_t m_init;
    double m_tol;

    coeff_type m_coeff;
};

}
}
}
}
#endif