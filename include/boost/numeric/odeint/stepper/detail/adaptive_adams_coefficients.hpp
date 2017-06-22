#ifndef ADAPTIVE_ADAMS_COEFFICIENTS_HPP_INCLUDED
#define ADAPTIVE_ADAMS_COEFFICIENTS_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>
#include <boost/numeric/odeint/stepper/detail/polynomial.hpp>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>
#include <boost/numeric/odeint/util/bind.hpp>

#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/operations_dispatcher.hpp>

#include <boost/numeric/odeint/util/unwrap_reference.hpp>

#include <iostream>
#include <boost/array.hpp>

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template<
size_t Steps,
class Deriv,
class Time,
class Algebra = typename algebra_dispatcher< Deriv >::algebra_type,
class Operations = typename operations_dispatcher< Deriv >::operations_type ,
class Resizer = initially_resizer
>
struct adaptive_adams_coefficients
{
public:
    static const size_t steps = Steps;

    typedef Deriv deriv_type;
    typedef Time time_type;
    typedef state_wrapper<deriv_type> wrapped_deriv_type;
    typedef detail::rotating_buffer<wrapped_deriv_type, steps+1> step_storage_type; // +1 for moulton
    typedef detail::rotating_buffer<time_type, steps+1> time_storage_type;

    typedef Algebra algebra_type;
    typedef Operations operations_type;

    typedef Resizer resizer_type;

    typedef adaptive_adams_coefficients<Steps, Deriv, Time, Algebra, Operations, Resizer> aac_type;
    typedef detail::Polynomial<steps+2, time_type> poly_type;

    adaptive_adams_coefficients(const algebra_type &algebra = algebra_type() )
    :poly(), m_effective_order(1), m_resizer(), m_algebra(algebra)
    {};

    void step(const deriv_type &deriv, const time_type &t)
    {
        m_resizer.adjust_size( deriv , detail::bind( &aac_type::template resize_tss_impl< deriv_type > , detail::ref( *this ) , detail::_1 ) );

        m_tts[0] = t;
        m_tss[0][0].m_v = deriv;

        for(size_t i=1; i<m_effective_order; ++i)
        {
            time_type dt = t - m_ts[i-1];
            m_algebra.for_each3(m_tss[i][0].m_v, m_tss[i-1][0].m_v, m_ss[i-1][0].m_v, typename Operations::template scale_sum2<double, double>(1/dt, -1/dt));
        }
    };
    void confirm()
    {
        for(size_t i=0; i<steps; ++i)
        {
            m_ss[i] = m_tss[i];
            m_tss[i].rotate();
        }

        m_ts = m_tts;
        m_tts.rotate();

        if(m_effective_order < steps+1)
        {
            ++m_effective_order;
        }
    };

    void reset()
    {
        poly.reset();
        m_effective_order = 1;
    };

    void pretty_print()
    {
        for(size_t k=0; k<2; ++k)
        {
            for(size_t i=0; i<m_effective_order; ++i)
            {
                for(size_t j=0; j<m_effective_order - i-1; ++j)
                    std::cout << m_ss[j][i].m_v[k] << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    };

    poly_type poly;
    time_storage_type m_c;

    boost::array<step_storage_type, steps+1> m_ss;
    time_storage_type m_ts;

    boost::array<step_storage_type, steps+1> m_tss;
    time_storage_type m_tts;

    size_t m_effective_order;

private:
    template< class StateType >
    bool resize_tss_impl( const StateType &x )
    {
        bool resized( false );
        for( size_t i=0 ; i<steps+1 ; ++i )
        {
            for(size_t j=0; j<steps+1; ++j)
                resized |= adjust_size_by_resizeability( m_tss[i][j] , x , typename is_resizeable<deriv_type>::type() );
            
            m_tts[i] = 0;
        }
        return resized;
    };

    resizer_type m_resizer;
    algebra_type m_algebra;
};

}
}
}
}

#endif
