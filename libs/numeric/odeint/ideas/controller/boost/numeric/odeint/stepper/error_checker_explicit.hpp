/*
 [auto_generated]
 boost/numeric/odeint/stepper/error_checker_explicit.hpp

 [begin_description]
 Error checker for the use in the generic_controlled_steppers in combination with explicit_error_steppers.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_ERROR_CHECKER_EXPLICIT_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_ERROR_CHECKER_EXPLICIT_HPP_INCLUDED


#include <boost/numeric/odeint/algebra/default_operations.hpp>


namespace boost {
namespace numeric {
namespace odeint {


/*
 * Error checker for controlled_error_stepper
 *
 * ToDo: remove, is redundant now
 */
template< class Value ,  class Algebra  , class Operations >
class error_checker_explicit
{
public:

    typedef Value value_type;
    typedef Algebra algebra_type;
    typedef Operations operations_type;


    error_checker_explicit(
            const value_type eps_abs = static_cast< value_type >( 1.0e-6 ) ,
            const value_type eps_rel = static_cast< value_type >( 1.0e-6 ) ,
            const value_type a_x = static_cast< value_type >( 1.0 ) ,
            const value_type a_dxdt = static_cast< value_type >( 1.0 ) )
    : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) , m_a_x( a_x ) , m_a_dxdt( a_dxdt )
    { }


    template< class State1 , class State2 , class Deriv , class Err , class Time >
    value_type error( const State1 &x_old , const State2 &x , const Deriv &dxdt_old , Err &x_err , const Time &dt )
    {
        // this overwrites x_err !
        algebra.for_each3( x_err , x_old , dxdt_old ,
                typename operations_type::template rel_error< value_type >( m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt * detail::get_value( dt ) ) );

        value_type res = algebra.reduce( x_err ,
                typename operations_type::template maximum< value_type >() , static_cast< value_type >( 0.0 ) );
        return res;
    }

private:

    value_type m_eps_abs;
    value_type m_eps_rel;
    value_type m_a_x;
    value_type m_a_dxdt;
    algebra_type algebra;
};





/*
 * Error checker for controlled_error_stepper
 *
 * ToDo: rename to error_checker_explicit_new
 */
template< class Value ,  class Algebra  , class Operations >
class error_checker_explicit_new
{
public:

    typedef Value value_type;
    typedef Algebra algebra_type;
    typedef Operations operations_type;


    error_checker_explicit_new(
            algebra_type &algebra ,
            const value_type eps_abs = static_cast< value_type >( 1.0e-6 ) ,
            const value_type eps_rel = static_cast< value_type >( 1.0e-6 ) ,
            const value_type a_x = static_cast< value_type >( 1.0 ) ,
            const value_type a_dxdt = static_cast< value_type >( 1.0 ) )
    : m_algebra( algebra ) , m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) , m_a_x( a_x ) , m_a_dxdt( a_dxdt )
    { }


    template< class State1 , class State2 , class Deriv , class Err , class Time >
    value_type error( const State1 &x_old , const State2 &x , const Deriv &dxdt_old , Err &x_err , const Time &dt )
    {
        // this overwrites x_err !
        m_algebra.for_each3( x_err , x_old , dxdt_old ,
                typename operations_type::template rel_error< value_type >( m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt * detail::get_value( dt ) ) );

        value_type res = m_algebra.reduce( x_err ,
                typename operations_type::template maximum< value_type >() , static_cast< value_type >( 0.0 ) );
        return res;
    }

private:

    algebra_type &m_algebra;
    value_type m_eps_abs;
    value_type m_eps_rel;
    value_type m_a_x;
    value_type m_a_dxdt;
};



}
}
}

#endif // BOOST_NUMERIC_ODEINT_STEPPER_ERROR_CHECKER_EXPLICIT_HPP_INCLUDED
