/*
 [auto_generated]
 boost/numeric/odeint/stepper/adams_bashforth.hpp

 [begin_description]
 Implementaton of the Adam-Bashforth method a multistep method used for the predictor step in the
 Adams-Bashforth-Moulton method.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_BASHFORTH_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_BASHFORTH_HPP_INCLUDED


#include <boost/ref.hpp>
#include <boost/bind.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

#include <boost/numeric/odeint/stepper/detail/adams_bashforth_coefficients.hpp>
#include <boost/numeric/odeint/stepper/detail/adams_bashforth_call_algebra.hpp>
#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>



namespace boost {
namespace numeric {
namespace odeint {


/*
 * Static explicit Adams-Bashforth multistep-solver without step size control and without dense output.
 */
template<
size_t Steps ,
class State ,
class Value = double ,
class Deriv = State ,
class Time = Value ,
class Algebra = range_algebra ,
class Operations = default_operations ,
class Resizer = initially_resizer
>
class adams_bashforth
{
private:

public :

    typedef State state_type;
    typedef state_wrapper< state_type > wrapped_state_type;
    typedef Value value_type;
    typedef Deriv deriv_type;
    typedef state_wrapper< deriv_type > wrapped_deriv_type;
    typedef Time time_type;
    typedef Algebra algebra_type;
    typedef Operations operations_type;
    typedef Resizer resizer_type;
    typedef stepper_tag stepper_category;

    typedef adams_bashforth< Steps , State , Value , Deriv , Time , Algebra , Operations , Resizer > stepper_type;

    static const size_t steps = Steps;

    typedef unsigned short order_type;
    static const order_type order_value = steps;

    typedef detail::rotating_buffer< wrapped_deriv_type , steps > step_storage_type;


    order_type order( void ) const { return order_value; }


    adams_bashforth( const algebra_type &algebra = algebra_type() ) : m_algebra( algebra )
    { }

    adams_bashforth( const adams_bashforth &stepper )
    : m_step_storage( stepper.m_step_storage ) , m_resizer( stepper.m_resizer ) , m_coefficients() , m_algebra( stepper.m_algebra )
    { }

    adams_bashforth& operator=( const adams_bashforth &stepper )
    {
        m_resizer = stepper.m_resizer;
        m_step_storage = stepper.m_step_storage;
        m_algebra = stepper.m_algebra;
        return *this;
    }


    /*
     * Version 1 : do_step( system , x , t , dt );
     *
     * solves the forwarding problem
     */
    template< class System , class StateInOut >
    void do_step( System system , StateInOut &x , const time_type &t , const time_type &dt )
    {
        do_step( system , x , t , x , dt );
    }

    template< class System , class StateInOut >
    void do_step( System system , const StateInOut &x , const time_type &t , const time_type &dt )
    {
        do_step( system , x , t , x , dt );
    }



    /*
     * Version 2 : do_step( system , in , t , out , dt );
     *
     * solves the forwarding problem
     */
    template< class System , class StateIn , class StateOut >
    void do_step( System system , const StateIn &in , const time_type &t , StateOut &out , const time_type &dt )
    {
        typename boost::unwrap_reference< System >::type &sys = system;
        m_resizer.adjust_size( in , boost::bind( &stepper_type::template resize<StateIn> , boost::ref( *this ) , _1 ) );
        m_step_storage.rotate();
        sys( in , m_step_storage[0].m_v , t );
        detail::adams_bashforth_call_algebra< steps , algebra_type , operations_type >()( m_algebra , in , out , m_step_storage , m_coefficients , dt );
    }

    template< class System , class StateIn , class StateOut >
    void do_step( System system , const StateIn &in , const time_type &t , const StateOut &out , const time_type &dt )
    {
        typename boost::unwrap_reference< System >::type &sys = system;
        m_resizer.adjust_size( in , boost::bind( &stepper_type::template resize<StateIn> , boost::ref( *this ) , _1 ) );
        m_step_storage.rotate();
        sys( in , m_step_storage[0].m_v , t );
        detail::adams_bashforth_call_algebra< steps , algebra_type , operations_type >()( m_algebra , in , out , m_step_storage , m_coefficients , dt );
    }




    //	/*
    //	 * Version 3 : do_step( system , x , dxdt , t , dt );
    //	 *
    //	 * solves the forwarding proble
    //	 *
    //	 * ToDo: Do we need this methods?
    //	 */
    //	template< class System , class StateInOut , class DerivIn >
    //	void do_step( System sys , StateInOut &x , const DerivIn &dxdt , const time_type &t , const time_type &dt )
    //	{
    //		do_step( sys , x , dxdt , t , x , dt );
    //	}
    //
    //	template< class System , class StateInOut , class DerivIn >
    //	void do_step( System sys , const StateInOut &x , const DerivIn &dxdt , const time_type &t , const time_type &dt )
    //	{
    //		do_step( sys , x , dxdt , t , x , dt );
    //	}
    //
    //
    //
    //	/*
    //	 * Version 4 : do_step( system , in , dxdt , t , out , dt )
    //	 *
    //	 * solves the forwarding problem
    //	 *
    // 	 * ToDo: Do we need this methods?
    //	 */
    //	template< class System , class StateIn , class DerivIn , class StateOut >
    //	void do_step( System system , const StateIn &in , const DerivIn &dxdt , const time_type &t , StateOut &out , const time_type &dt )
    //	{
    //		m_step_storage.rotate();
    //		boost::numeric::odeint::copy( dxdt , m_step_storage[0] );
    //		do_step_impl( in , t , out , dt );
    //	}
    //
    //	template< class System , class StateIn , class DerivIn , class StateOut >
    //	void do_step( System system , const StateIn &in , const DerivIn &dxdt , const time_type &t , const StateOut &out , const time_type &dt )
    //	{
    //		m_step_storage.rotate();
    //		boost::numeric::odeint::copy( dxdt , m_step_storage[0] );
    //		do_step_impl( in , t , out , dt );
    //	}










    template< class StateIn >
    bool resize( const StateIn &x )
    {
        bool resized( false );
        for( size_t i=0 ; i<steps ; ++i )
        {
            resized |= adjust_size_by_resizeability( m_step_storage[i] , x , typename wrapped_deriv_type::is_resizeable() );
        }
        return resized;
    }


    template< class StateType >
    void adjust_size( const StateType &x )
    {
        resize( x );
    }

    algebra_type& algebra()
    {   return m_algebra; }

    const algebra_type& algebra() const
    {   return m_algebra; }

    const step_storage_type& step_storage( void ) const
    {
        return m_step_storage;
    }

    step_storage_type& step_storage( void )
    {
        return m_step_storage;
    }

    template< class ExplicitStepper , class System , class StateIn >
    void initialize( ExplicitStepper explicit_stepper , System system , StateIn &x , time_type &t , const time_type &dt )
    {
        typename boost::unwrap_reference< ExplicitStepper >::type &stepper = explicit_stepper;
        typename boost::unwrap_reference< System >::type &sys = system;

        m_resizer.adjust_size( x , boost::bind( &stepper_type::template resize<StateIn> , boost::ref( *this ) , _1 ) );

        for( size_t i=0 ; i<steps-1 ; ++i )
        {
            if( i != 0 ) m_step_storage.rotate();
            sys( x , m_step_storage[0].m_v , t );
            stepper.do_step( system , x , m_step_storage[0].m_v , t , dt );
            t += dt;
        }
    }

    template< class System , class StateIn >
    void initialize( System system , StateIn &x , time_type &t , const time_type &dt )
    {
        /**
         * ToDo: rk4 here has its own algebra instance! good? bad?
         */
        runge_kutta4< state_type , value_type , deriv_type , time_type , algebra_type , operations_type , resizer_type > rk4;
        initialize( boost::ref( rk4 ) , system , x , t , dt );
    }


private:

    step_storage_type m_step_storage;
    resizer_type m_resizer;
    const detail::adams_bashforth_coefficients< value_type , steps > m_coefficients;

protected:
    algebra_type m_algebra;
};




} // odeint
} // numeric
} // boost



#endif // BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_BASHFORTH_HPP_INCLUDED
