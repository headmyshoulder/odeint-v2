#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_ADAPTIVE_ADAMS_COEFFICIENTS_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_ADAPTIVE_ADAMS_COEFFICIENTS_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/numeric/odeint/util/bind.hpp>

#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/operations_dispatcher.hpp>

#include <boost/array.hpp>

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template<
size_t Steps,
class Deriv,
class Value = double,
class Time = double,
class Algebra = typename algebra_dispatcher< Deriv >::algebra_type,
class Operations = typename operations_dispatcher< Deriv >::operations_type,
class Resizer = initially_resizer
>
class adaptive_adams_coefficients
{
public:
    static const size_t steps = Steps;

    typedef unsigned short order_type;
    static const order_type order_value = steps;

    typedef Value value_type;
    typedef Deriv deriv_type;
    typedef Time time_type;

    typedef state_wrapper< deriv_type > wrapped_deriv_type;
    typedef rotating_buffer< wrapped_deriv_type , steps+1 > step_storage_type; // +1 for moulton
    typedef rotating_buffer< time_type , steps+1 > time_storage_type;

    typedef Algebra algebra_type;
    typedef Operations operations_type;
    typedef Resizer resizer_type;

    typedef adaptive_adams_coefficients< Steps , Deriv , Value , Time , Algebra , Operations , Resizer > aac_type;

    adaptive_adams_coefficients( const algebra_type &algebra = algebra_type())
    :m_eo(1), m_steps_init(1), beta(), phi(), m_ns(0), m_time_storage(),
    m_algebra(algebra),
    m_phi_resizer()
    {
        for (size_t i=0; i<order_value+2; ++i)
        {
            c[0][i] = 1.0/(i+1);
        }

        g[0] = c[0][0];

        beta[0][0] = 1;
        beta[1][0] = 1;
    };

    void predict(time_type t, time_type dt)
    {
        using std::abs;

        m_time_storage[0] = t;

        if (abs(m_time_storage[0] - m_time_storage[1] - dt) > 1e-16 || m_eo >= m_ns)
        {
            m_ns = 0;
        }
        else if (m_ns < order_value + 2)
        {
            m_ns++;
        }

        for(size_t i=1+m_ns; i<m_eo+1 && i<m_steps_init; ++i)
        {
            beta[0][i] = beta[0][i-1]*(m_time_storage[0] + dt -
                m_time_storage[i-1])/(m_time_storage[0] - m_time_storage[i]);
        }

        for(size_t i=1+m_ns; i<m_eo+2 && i<m_steps_init+1; ++i)
        {
            for(size_t j=0; j<m_eo+1; ++j)
            {
                c[i][j] = c[i-1][j] - c[i-1][j+1]*dt/(m_time_storage[0] + dt - m_time_storage[i-1]);
            }

            g[i] = c[i][0];
        }
    };

    void do_step(const deriv_type &dxdt, const int o = 0)
    {
        m_phi_resizer.adjust_size( dxdt , detail::bind( &aac_type::template resize_phi_impl< deriv_type > , detail::ref( *this ) , detail::_1 ) );

        phi[o][0].m_v = dxdt;

        for(size_t i=1; i<m_eo+2 && i<m_steps_init+1; ++i)
        {
            this->m_algebra.for_each3(phi[o][i].m_v, phi[o][i-1].m_v, phi[o+1][i-1].m_v,
                typename Operations::template scale_sum2<double, double>(1.0, -beta[o][i-1]));
        }   
    };

    void confirm()
    {
        beta.rotate();
        phi.rotate();
        m_time_storage.rotate();

        if(m_steps_init < order_value+1)
        {
            ++m_steps_init;
        }
    };

    void reset() { m_eo = 1; };

    size_t m_eo;
    size_t m_steps_init;

    rotating_buffer<boost::array<value_type, order_value+1>, 2> beta; // beta[0] = beta(n)
    rotating_buffer<boost::array<wrapped_deriv_type, order_value+2>, 3> phi; // phi[0] = phi(n+1)
    boost::array<value_type, order_value + 2> g;

private:
    template< class StateType >
    bool resize_phi_impl( const StateType &x )
    {

        bool resized( false );

        for(size_t i=0; i<(order_value + 2); ++i)
        {
            resized |= adjust_size_by_resizeability( phi[0][i], x, typename is_resizeable<deriv_type>::type() );
            resized |= adjust_size_by_resizeability( phi[1][i], x, typename is_resizeable<deriv_type>::type() );
            resized |= adjust_size_by_resizeability( phi[2][i], x, typename is_resizeable<deriv_type>::type() );
        }
        return resized;
    };

    size_t m_ns;

    time_storage_type m_time_storage;
    boost::array<boost::array<value_type, order_value + 2>, order_value + 2> c;

    algebra_type m_algebra;

    resizer_type m_phi_resizer;
};

} // detail
} // odeint
} // numeric
} // boost

#endif