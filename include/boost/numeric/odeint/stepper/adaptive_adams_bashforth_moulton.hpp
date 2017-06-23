#ifndef STEPPER_ADAPTIVE_ADAMS_BASHFORTH_MOULTON_HPP_INCLUDED
#define STEPPER_ADAPTIVE_ADAMS_BASHFORTH_MOULTON_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/base/algebra_stepper_base.hpp>
#include <boost/numeric/odeint/stepper/adaptive_adams_bashforth.hpp>
#include <boost/numeric/odeint/stepper/adaptive_adams_moulton.hpp>
#include <boost/numeric/odeint/stepper/detail/adaptive_adams_coefficients.hpp>

#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/operations_dispatcher.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/util/bind.hpp>

namespace boost {
namespace numeric {
namespace odeint {



template<
size_t Steps,
class State,
class Value = double,
class Deriv = State,
class Time = Value,
class Algebra = typename algebra_dispatcher< State >::algebra_type,
class Operations = typename operations_dispatcher< State >::operations_type ,
class Resizer = initially_resizer
>
class adaptive_adams_bashforth_moulton: public algebra_stepper_base< Algebra , Operations >
{
public:
    static const size_t steps = Steps;
    typedef unsigned short order_type;
    static const order_type order_value = steps;

    typedef State state_type;
    typedef Value value_type;
    typedef Deriv deriv_type;
    typedef Time time_type;
    typedef Resizer resizer_type;

    typedef Algebra algebra_type;
    typedef Operations operations_type;
    typedef algebra_stepper_base< Algebra , Operations > algebra_stepper_base_type;

    typedef state_wrapper<state_type> wrapped_state_type;
    typedef state_wrapper<deriv_type> wrapped_deriv_type;

    typedef detail::adaptive_adams_coefficients<order_value, deriv_type, time_type, algebra_type, operations_type> coeff_type;

    typedef detail::rotating_buffer<state_type, steps> error_storage_type;

    typedef adaptive_adams_bashforth<order_value, state_type, value_type, state_type, value_type, algebra_type, operations_type> aab_type;
    typedef adaptive_adams_moulton<order_value, state_type, value_type, state_type, value_type, algebra_type, operations_type> aam_type;
    typedef adaptive_adams_bashforth_moulton< Steps , State , Value , Deriv , Time, Algebra, Operations, Resizer > stepper_type;

    typedef error_stepper_tag stepper_category;

    adaptive_adams_bashforth_moulton(const algebra_type &algebra = algebra_type())
    :algebra_stepper_base_type( algebra ), m_aab(algebra), m_aam(algebra),
    m_dxdt_resizer(), m_xerr_resizer(), m_xnew_resizer(),
    m_coeff()
    {};

    order_type order()
    {
        return order_value;
    }

    order_type stepper_order()
    {
        return order_value + 1;
    }

    order_type error_order()
    {
        return order_value;
    }

    template<class System>
    void do_step(System system, state_type & inOut, time_type t, time_type &dt)
    {
        m_xerr_resizer.adjust_size( inOut , detail::bind( &stepper_type::template resize_xerr_impl< state_type > , detail::ref( *this ) , detail::_1 ) );

        do_step(system, inOut, t, dt, m_xerr.m_v);
    };

    template<class System>
    void do_step(System system, const state_type & in, time_type t, state_type & out, time_type &dt)
    {
        m_xerr_resizer.adjust_size( in , detail::bind( &stepper_type::template resize_xerr_impl< state_type > , detail::ref( *this ) , detail::_1 ) );

        do_step(system, in, t, out, dt, m_xerr.m_v);
    };

    template<class System>
    void do_step(System system, state_type & inOut, time_type t, time_type &dt, state_type &xerr)
    {
        m_xnew_resizer.adjust_size( inOut , detail::bind( &stepper_type::template resize_xnew_impl< state_type > , detail::ref( *this ) , detail::_1 ) );

        do_step(system, inOut, t, m_xnew.m_v, dt, xerr);
        boost::numeric::odeint::copy( m_xnew.m_v , inOut);
    };

    template<class System>
    void do_step(System system, const state_type & in, time_type &t, state_type & out, time_type dt, state_type &xerr)
    {
        do_step_impl(system, m_coeff, in, t, out, dt, xerr);

        system(out, m_dxdt.m_v, t+dt);
        m_coeff.step(m_dxdt.m_v, t+dt);
        m_coeff.confirm();
    };

    template<class System>
    void do_step_impl(System system, coeff_type coeff, const state_type & in, time_type &t, state_type & out, time_type &dt, state_type &xerr)
    {
        m_dxdt_resizer.adjust_size( in , detail::bind( &stepper_type::template resize_dxdt_impl< state_type > , detail::ref( *this ) , detail::_1 ) );
        
        if(coeff.m_effective_order == 1)
        {
            system(in, m_dxdt.m_v, t);

            coeff.step(m_dxdt.m_v, t);
            coeff.confirm();
        }
        // predict
        m_aab.do_step_impl(coeff, in, t, out, dt);

        // evaluate
        system(out, m_dxdt.m_v, t + dt);
        coeff.step(m_dxdt.m_v, t + dt);

        boost::numeric::odeint::copy( out, xerr);

        // correct
        m_aam.do_step(coeff, in, t, out, dt);

        // Error calculation
        m_aab.algebra().for_each3(xerr, xerr, out, typename Operations::template scale_sum2<double, double>(-1.0, 1.0));
    };

    template<class System>
    void initialize(System system, state_type &inOut, time_type &t, time_type dt)
    {
        m_coeff.reset();

        for(size_t i=0; i<steps+1; ++i)
        {
            do_step(system, inOut, t, dt);
            t += dt;
        }
    }

    const coeff_type& coeff() const
    {
        return m_coeff;
    };

    coeff_type& coeff()
    {
        return m_coeff;
    };

    wrapped_deriv_type m_dxdt;

private:
    template< class StateType >
    bool resize_dxdt_impl( const StateType &x )
    {
        return adjust_size_by_resizeability( m_dxdt, x, typename is_resizeable<deriv_type>::type() );
    };
    template< class StateType >
    bool resize_xerr_impl( const StateType &x )
    {
        return adjust_size_by_resizeability( m_xerr, x, typename is_resizeable<deriv_type>::type() );
    };
    template< class StateType >
    bool resize_xnew_impl( const StateType &x )
    {
        return adjust_size_by_resizeability( m_xnew, x, typename is_resizeable<state_type>::type() );
    };

    aab_type m_aab;
    aam_type m_aam;

    resizer_type m_dxdt_resizer;
    resizer_type m_xerr_resizer;
    resizer_type m_xnew_resizer;

    wrapped_state_type m_xnew;
    wrapped_state_type m_xerr;

    coeff_type m_coeff;
};

}}}

#endif
