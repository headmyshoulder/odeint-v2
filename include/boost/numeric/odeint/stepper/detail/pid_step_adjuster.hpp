#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_PID_STEP_ADJUSTER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_PID_STEP_ADJUSTER_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>
#include <boost/numeric/odeint/stepper/detail/pid_step_adjuster_coefficients.hpp>

#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>

#include <math.h>

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template<
class Value = double,
class Time = double
>
struct pid_op
{
public:
    typedef Value value_type;
    typedef Time time_type;

    const double beta1;
    const double beta2;
    const double beta3;
    const double alpha1;
    const double alpha2;

    const time_type dt1;
    const time_type dt2;
    const time_type dt3;

    const double m_tol;
    const size_t m_steps;

    pid_op(const size_t steps, const double tol, const double _dt1, const double _dt2, const double _dt3,
        const double b1 = 1, const double b2 = 0, const double b3 = 0, const double a1 = 0, const double a2 = 0)
    :beta1(b1), beta2(b2), beta3(b3), alpha1(a1), alpha2(a2),
    dt1(_dt1), dt2(_dt2), dt3(_dt3),
    m_tol(tol), m_steps(steps)
    {};

    template<class T1, class T2, class T3, class T4>
    void operator()(T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4)
    {
        t1 = adapted_pow(fabs(m_tol/t2), beta1/(m_steps + 1)) *
            adapted_pow(fabs(m_tol/t3), beta2/(m_steps + 1)) *
            adapted_pow(fabs(m_tol/t4), beta3/(m_steps + 1)) *
            adapted_pow(fabs(dt1/dt2), -alpha1/(m_steps + 1))*
            adapted_pow(fabs(dt2/dt3), -alpha2/(m_steps + 1));

        t1 = 1/t1;
    };

    template<class T1, class T2>
    void operator()(T1 &t1, const T2 &t2)
    {
        t1 = adapted_pow(fabs(m_tol/t2), beta1/(m_steps + 1));

        t1 = 1/t1;
    };

private:
    template<class T>
    inline value_type adapted_pow(T base, double exp)
    {
        if(exp == 0)
        {
            return 1;
        }
        else if (exp > 0)
        {
            return pow(base, exp);
        }
        else
        {
            return 1/pow(base, -exp);
        }
    };
};

template<
class State,
class Value = double,
class Time = double,
class Algebra = typename algebra_dispatcher< State >::algebra_type,
size_t Type = BASIC
>
struct pid_step_adjuster
{
public:
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

    time_type adjust_stepsize(const size_t steps, time_type dt, state_type &err)
    {
        m_error_storage[0] = err;
        m_dt_storage[0] = dt;

        if(m_init >= 2)
        {
            m_algebra.for_each4(err, m_error_storage[0], m_error_storage[1], m_error_storage[2],
                pid_op<>(steps, m_tol, m_dt_storage[0], m_dt_storage[1], m_dt_storage[2],
                    m_coeff[0], m_coeff[1], m_coeff[2], m_coeff[3], m_coeff[4]));
        }
        else
        {
            m_algebra.for_each2(err, m_error_storage[0],
                pid_op<>(steps, m_tol, m_dt_storage[0], m_dt_storage[1], m_dt_storage[2], 0.7));
        }

        value_type ratio = 1 / m_algebra.norm_inf(err);

        value_type kappa = 1.0;
        ratio = 1.0 + kappa*atan((ratio - 1) / kappa);

        if(ratio*dt >= m_dtmax)
        {
            ratio = m_dtmax / dt;
        }

        if(ratio >= threshold())
        {
            m_error_storage.rotate();
            m_dt_storage.rotate();

            ++m_init;
        }
        else
        {
            m_init = 0;
        }

        return dt * static_cast<time_type>(ratio);
    };

private:
    algebra_type m_algebra;

    time_type m_dtmax;
    error_storage_type m_error_storage;
    time_storage_type m_dt_storage;

    size_t m_init;
    double m_tol;

    coeff_type m_coeff;
};

} // detail
} // odeint
} // numeric
} // boost

#endif