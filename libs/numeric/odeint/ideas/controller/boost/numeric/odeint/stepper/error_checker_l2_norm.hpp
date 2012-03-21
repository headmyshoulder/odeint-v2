/*
 [auto_generated]
 boost/numeric/odeint/stepper/default_error_checker_l2_norm.hpp

 [begin_description]
 Default error checker for the use in the generic_controlled_steppers. Works with all error_steppers, explicit_error_steppers and
 explicit_error_stepper_fsal;
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_ERROR_CHECKER_L2_NORM_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_ERROR_CHECKER_L2_NORM_HPP_INCLUDED


#include <boost/numeric/odeint/util/ref_or_value_holder.hpp>
#include <boost/numeric/odeint/util/number_of_elements.hpp>

namespace boost {
namespace numeric {
namespace odeint {



/*
 * Calculate the error in a generic way:
 *
 * for steppers: err = ||( | err[i] | / ( eps_abs + eps_rel * | x[i] | ) )||
 * for explicit steppers: err = ||( | err[i] | ( eps_abs + eps_rel * ( a_x * | x[i] | + a_dxdt * | dxdt[i] * dt ] ) ) )||
 * for explicit fsal steppers: err = ||( | err[i] | ( eps_abs + eps_rel * ( a_x * | x[i] | + a_dxdt * | dxdt[i] * dt ] ) ) )||
 */
template< class Value ,  class Algebra  , class Operations , bool HoldsAlgebra = true >
class error_checker_l2_norm
{
public:

    typedef Value value_type;
    typedef Algebra algebra_type;
    typedef Operations operations_type;
    const static bool holds_algebra = HoldsAlgebra;
    const static bool ref_algebra = !holds_algebra;
    typedef ref_or_value_holder< algebra_type , ref_algebra > algebra_holder_type;

    error_checker_l2_norm(
            typename algebra_holder_type::constructor_type algebra ,
            const value_type eps_abs = static_cast< value_type >( 1.0e-6 ) ,
            const value_type eps_rel = static_cast< value_type >( 1.0e-6 ) ,
            const value_type a_x = static_cast< value_type >( 1.0 ) ,
            const value_type a_dxdt = static_cast< value_type >( 1.0 ))
    : m_algebra( algebra ) ,
      m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) , m_a_x( a_x ) , m_a_dxdt( a_dxdt )
    { }


    // overload for steppers, x, x_old and x_err are available
    template< class State1 , class State2 , class Err , class Time >
    value_type error( const State1 &x_old , const State2 &x , const Err &x_err , const Time &dt )
    {
        value_type result = m_algebra.get().reduce3( x_old , x , x_err ,
                typename operations_type::template rel_error_l2_2< value_type >( m_eps_abs , m_eps_rel ) ,
                static_cast< value_type >( 0.0 ) );
        result /= value_type( boost::numeric::odeint::size( x_old ) );
        return result;
    }


    // overload for explicit steppers, x, x_old, dxdt_old and x_err are available
    template< class State1 , class State2 , class Deriv , class Err , class Time >
    value_type error( const State1 &x_old , const State2 &x , const Deriv &dxdt_old , Err &x_err , const Time &dt )
    {
        value_type result = m_algebra.get().reduce4( x_old , x , dxdt_old , x_err ,
                typename operations_type::template rel_error_l2_2< value_type >( m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt * detail::get_value( dt ) ) ,
                static_cast< value_type >( 0.0 ) );
        result /= value_type( boost::numeric::odeint::size( x_old ) );
        return result;
    }


    // overload for explicit fsal steppers, x, x_old, dxdt, dxdt_old and x_err are available
    template< class StateOld , class State , class DerivOld , class Deriv , class Err , class Time >
    value_type error( const StateOld &x_old , const State &x , const DerivOld &dxdt_old , const Deriv &dxdt , const Err &x_err , const Time &dt )
    {
        value_type result = m_algebra.get().reduce4( x_old , x , dxdt_old , x_err ,
                typename operations_type::template rel_error_l2_2< value_type >( m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt * detail::get_value( dt ) ) ,
                static_cast< value_type >( 0.0) );
        result /= value_type( boost::numeric::odeint::size( x_old ) );
        return result;
    }



private:

    algebra_holder_type m_algebra;
    value_type m_eps_abs;
    value_type m_eps_rel;
    value_type m_a_x;
    value_type m_a_dxdt;
};


}
}
}

#endif // BOOST_NUMERIC_ODEINT_STEPPER_ERROR_CHECKER_EXPLICIT_L2_NORM_HPP_INCLUDED
