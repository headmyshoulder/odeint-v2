/*
 [auto_generated]
 boost/numeric/odeint/stepper/generic_controlled_stepper.hpp

 [begin_description]
 Implementation of the generic_controlled_stepper for the error_stepper_tag. Specialized versions
 for explicit_error_stepper and explicit_error_stepper_tag also exist.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_GENERIC_CONTROLLED_STEPPER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_GENERIC_CONTROLLED_STEPPER_HPP_INCLUDED


#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>

#include <boost/numeric/odeint/stepper/generic_controlled_stepper_definition.hpp>

#include <boost/bind.hpp>



namespace boost {
namespace numeric {
namespace odeint {




template<
    class ErrorStepper ,
    class ErrorChecker ,
    class Controller ,
    class Resizer >
class generic_controlled_stepper< ErrorStepper , ErrorChecker , Controller , Resizer , error_stepper_tag >
{
public:

    typedef ErrorStepper stepper_type;
    typedef ErrorChecker error_checker_type;
    typedef Controller controller_type;
    typedef Resizer resizer_type;
    typedef typename stepper_type::state_type state_type;
    typedef typename stepper_type::value_type value_type;
    typedef typename stepper_type::time_type time_type;
    typedef typename stepper_type::order_type order_type;
    typedef explicit_controlled_stepper_tag stepper_category;
    typedef state_wrapper< state_type > wrapped_state_type;



    generic_controlled_stepper(
            const stepper_type &stepper = stepper_type() ,
            const error_checker_type &error_checker = error_checker_type() ,
            const controller_type &controller = controller_type() )
    : m_stepper( stepper ) , m_error_checker( error_checker ) , m_controller( controller ) ,
      m_xerr_resizer() , m_xnew_resizer() ,
      m_xerr() , m_xnew() { }

    /*
    * Version 1 : try_step( sys , x , t , dt )
    *
    * The overloads are needed to solve the forwarding problem
    */
    template< class System , class StateInOut >
    controlled_step_result try_step( System system , StateInOut &x , time_type &t , time_type &dt )
    {
        m_xnew_resizer.adjust_size( x , boost::bind( &generic_controlled_stepper::template resize_m_xnew_impl< state_type > , boost::ref( *this ) , _1 ) );
        controlled_step_result result = try_step( system , x , t , m_xnew.m_v , dt );
        if( result == success )
        {
            boost::numeric::odeint::copy( m_xnew.m_v , x );
        }
        return result;
    }

    template< class System , class StateInOut >
    controlled_step_result try_step( System system , const StateInOut &x , time_type &t , time_type &dt )
    {
        m_xnew_resizer.adjust_size( x , boost::bind( &generic_controlled_stepper::template resize_m_xnew_impl< state_type > , boost::ref( *this ) , _1 ) );
        controlled_step_result result = try_step( system , x , t , m_xnew.m_v , dt );
        if( result == success )
        {
            boost::numeric::odeint::copy( m_xnew.m_v , x );
        }
        return result;
    }



    /*
     * Version 2 : try_step( sys , in , t , out , dt )
     *
     * this version does not solve the forwarding problem, boost.range can not be used
     */
    template< class System , class StateIn , class StateOut >
    controlled_step_result try_step( System system , const StateIn &in , time_type &t , StateOut &out , time_type &dt )
    {
//        return try_step( system , in , m_dxdt.m_v , t , out , dt );
        // ToDo : implement

        m_xerr_resizer.adjust_size( in , boost::bind( &generic_controlled_stepper::template resize_m_xerr_impl< StateIn > , boost::ref( *this ) , _1 ) );
        m_stepper.do_step( system , in , t , out , dt , m_xerr.m_v );
        value_type error = m_error_checker.error( in , out , m_xerr.m_v , dt );
        return m_controller( error , t , dt , m_stepper.stepper_order , m_stepper.error_order );

    }








    template< class StateType >
    void adjust_size( const StateType &x )
    {
        resize_m_xerr_impl( x );
        resize_m_xnew_impl( x );
        m_stepper.adjust_size( x );
    }


    stepper_type& stepper( void ) { return m_stepper; }
    const stepper_type stepper( void ) const { return m_stepper; }


private:


    template< class StateIn >
    bool resize_m_xerr_impl( const StateIn &x )
    {
        return adjust_size_by_resizeability( m_xerr , x , typename is_resizeable< wrapped_state_type >::type() );
    }

    template< class StateIn >
    bool resize_m_xnew_impl( const StateIn &x )
    {
        return adjust_size_by_resizeability( m_xnew , x , typename is_resizeable< wrapped_state_type >::type() );
    }




    stepper_type m_stepper;
    error_checker_type m_error_checker;
    controller_type m_controller;

    resizer_type m_xerr_resizer;
    resizer_type m_xnew_resizer;

    wrapped_state_type m_xerr;
    wrapped_state_type m_xnew;

};




}
}
}



#endif // BOOST_NUMERIC_ODEINT_STEPPER_GENERIC_CONTROLLED_STEPPER_HPP_INCLUDED
