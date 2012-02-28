/*
 [auto_generated]
 boost/numeric/odeint/stepper/default_error_checker.hpp

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


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DEFAULT_ERROR_CHECKER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_DEFAULT_ERROR_CHECKER_HPP_INCLUDED


#include <boost/numeric/odeint/algebra/default_operations.hpp>


namespace boost {
namespace numeric {
namespace odeint {



/*
 * Calculate the error in a generic way:
 *
 * for steppers: err = max_i( | err[i] | / ( eps_abs + eps_rel * | x[i] | ) )
 * for explicit steppers: err = max_i( | err[i] | ( eps_abs + eps_rel * ( a_x * | x[i] | + a_dxdt * | dxdt[i] * dt ] ) ) )
 * for explicit fsal steppers: err = max_i( | err[i] | ( eps_abs + eps_rel * ( a_x * | x[i] | + a_dxdt * | dxdt[i] * dt ] ) ) )
 */
template< class Value ,  class Algebra  , class Operations >
class default_error_checker
{
public:

    typedef Value value_type;
    typedef Algebra algebra_type;
    typedef Operations operations_type;


    default_error_checker( void )
    { }

    // overload for steppers, x, x_old and x_err are available
    template< class State1 , class State2 , class Err , class Time >
    value_type error( const State1 &x_old , const State2 &x , const Err &x_err , const Time &dt )
    {
        value_type result = reduce3( x_old , x , x_err ,
                typename operations_type::template rel_error_max< value_type >( m_eps_abs , m_eps_rel ) );
        return result;
    }

    // overload for explicit steppers, x, x_old, dxdt_old and x_err are available
    template< class State1 , class State2 , class Deriv , class Err , class Time >
    value_type error( const State1 &x_old , const State2 &x , const Deriv &dxdt_old , Err &x_err , const Time &dt )
    {
        value_type result = reduce4( x_old , x , dxdt_old , x_err ,
                typename operations_type::template rel_error_max2< value_type >( m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt ) );
        return result;
    }

    // overload for explicit fsal steppers, x, x_old, dxdt, dxdt_old and x_err are available
    template< class StateOld , class State , class DerivOld , class Deriv , class Err , class Time >
    value_type error( const StateOld &x_old , const State &x , const DerivOld &dxdt_old , const Deriv &dxdt , const Err &x_err , const Time &dt )
    {
        value_type result = reduce4( x_old , x , dxdt_old , x_err ,
                typename operations_type::template rel_error_max2< value_type >( m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt ) );
        return result;
    }



private:

    algebra_type algebra;
};


}
}
}

#endif // BOOST_NUMERIC_ODEINT_STEPPER_ERROR_CHECKER_EXPLICIT_HPP_INCLUDED
