/*
 boost header: numeric/odeint/explicit_stepper_and_error_stepper_fsal_base.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_EXPLICIT_STEPPER_AND_ERROR_STEPPER_FSAL_BASE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_EXPLICIT_STEPPER_AND_ERROR_STEPPER_FSAL_BASE_HPP_INCLUDED

#include <boost/ref.hpp>

#include <boost/numeric/odeint/util/size_adjuster.hpp>
#include <boost/numeric/odeint/util/construct.hpp>
#include <boost/numeric/odeint/util/destruct.hpp>
#include <boost/numeric/odeint/util/copy.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

namespace boost {
namespace numeric {
namespace odeint {

/*
 * base class for explicit stepper and error steppers with the fsal property
 * models the stepper AND the error stepper fsal concept
 */
template<
	class Stepper ,
	unsigned short Order ,
	unsigned short StepperOrder ,
	unsigned short ErrorOrder ,
	class State ,
	class Value ,
	class Deriv ,
	class Time ,
 	class Algebra ,
	class Operations ,
	class AdjustSizePolicy
>
class explicit_stepper_and_error_stepper_fsal_base
{
public:

	typedef State state_type;
	typedef Value value_type;
	typedef Deriv deriv_type;
	typedef Time time_type;
	typedef Algebra algebra_type;
	typedef Operations operations_type;
	typedef AdjustSizePolicy adjust_size_policy;
	typedef Stepper stepper_type;
	typedef explicit_error_stepper_fsal_tag stepper_category;

	typedef unsigned short order_type;
	static const order_type order_value = Order;
	static const order_type stepper_order_value = StepperOrder;
	static const order_type error_order_value = ErrorOrder;




    order_type order( void ) const
    {
    	return order_value;
    }

    order_type stepper_order( void ) const
    {
    	return stepper_order_value;
    }

    order_type error_order( void ) const
    {
    	return error_order_value;
    }




    explicit_stepper_and_error_stepper_fsal_base( void ) : m_size_adjuster() , m_dxdt() , m_first_call( true )
	{
		boost::numeric::odeint::construct( m_dxdt );
		m_size_adjuster.register_state( 0 , m_dxdt );
	}

	~explicit_stepper_and_error_stepper_fsal_base( void )
	{
		boost::numeric::odeint::destruct( m_dxdt );
	}

    explicit_stepper_and_error_stepper_fsal_base( const explicit_stepper_and_error_stepper_fsal_base &b )
    : m_size_adjuster() , m_dxdt() , m_first_call( true )
	{
		boost::numeric::odeint::construct( m_dxdt );
		m_size_adjuster.register_state( 0 , m_dxdt );
		boost::numeric::odeint::copy( b.m_dxdt , m_dxdt );
	}

    explicit_stepper_and_error_stepper_fsal_base& operator=( const explicit_stepper_and_error_stepper_fsal_base &b )
    {
    	boost::numeric::odeint::copy( b.m_dxdt , m_dxdt );
    	m_first_call( true );
		return *this;
    }



    /*
     * version 1 : do_step( sys , x , t , dt )
     *
     * the two overloads are needed in order to solve the forwarding problem
     */
	template< class System , class StateInOut >
	void do_step( System system , StateInOut &x , const time_type &t , const time_type &dt )
	{
		do_step_v1( system , x , t , dt );
	}

	template< class System , class StateInOut >
	void do_step( System system , const StateInOut &x , const time_type &t , const time_type &dt )
	{
		do_step_v1( system , x , t , dt );
	}


    /*
     * version 2 : do_step( sys , x , dxdt , t , dt )
     *
     * this version does not solve the forwarding problem, boost.range can not be used
     */
	template< class System , class StateInOut , class DerivInOut >
	void do_step( System system , StateInOut &x , DerivInOut &dxdt , const time_type &t , const time_type &dt )
	{
		m_first_call = true;
		this->stepper().do_step_impl( system , x , dxdt , t , x , dxdt , dt );
	}


    /*
     * version 3 : do_step( sys , in , t , out , dt )
     *
     * this version does not solve the forwarding problem, boost.range can not be used
     */
	template< class System , class StateIn , class StateOut >
	void do_step( System system , const StateIn &in , const time_type &t , StateOut &out , const time_type &dt )
	{
		if( m_size_adjuster.adjust_size_by_policy( in , adjust_size_policy() ) || m_first_call )
		{
			typename boost::unwrap_reference< System >::type &sys = system;
			sys( in , m_dxdt ,t );
			m_first_call = false;
		}
		this->stepper().do_step_impl( system , in , m_dxdt , t , out , m_dxdt , dt );
	}


    /*
     * version 4 : do_step( sys , in , dxdt_in , t , out , dxdt_out , dt )
     *
     * this version does not solve the forwarding problem, boost.range can not be used
     */
	template< class System , class StateIn , class DerivIn , class StateOut , class DerivOut >
	void do_step( System system , const StateIn &in , const DerivIn &dxdt_in , const time_type &t ,
			StateOut &out , DerivOut &dxdt_out , const time_type &dt )
	{
		m_first_call = true;
		this->stepper().do_step_impl( system , in , dxdt_in , t , out , dxdt_out , dt );
	}





    /*
     * version 5 : do_step( sys , x , t , dt , xerr )
     *
     * the two overloads are needed in order to solve the forwarding problem
     */
	template< class System , class StateInOut , class Err >
	void do_step( System system , StateInOut &x , const time_type &t , const time_type &dt , Err &xerr )
	{
		do_step_v5( system , x , t , dt , xerr );
	}

	template< class System , class StateInOut , class Err >
	void do_step( System system , const StateInOut &x , const time_type &t , const time_type &dt , Err &xerr )
	{
		do_step_v5( system , x , t , dt , xerr );
	}


    /*
     * version 6 : do_step( sys , x , dxdt , t , dt , xerr )
     *
     * this version does not solve the forwarding problem, boost.range can not be used
     */
	template< class System , class StateInOut , class DerivInOut , class Err >
	void do_step( System system , StateInOut &x , DerivInOut &dxdt , const time_type &t , const time_type &dt , Err &xerr )
	{
		m_first_call = true;
		this->stepper().do_step_impl( system , x , dxdt , t , x , dxdt , dt , xerr );
	}


    /*
     * version 7 : do_step( sys , in , t , out , dt , xerr )
     *
     * this version does not solve the forwarding problem, boost.range can not be used
     */
	template< class System , class StateIn , class StateOut , class Err >
	void do_step( System system , const StateIn &in , const time_type &t , StateOut &out , const time_type &dt , Err &xerr )
	{
	    if( m_size_adjuster.adjust_size_by_policy( in , adjust_size_policy() ) || m_first_call )
	    {
	    	typename boost::unwrap_reference< System >::type &sys = system;
	        sys( in , m_dxdt ,t );
	        m_first_call = false;
	    }
		this->stepper().do_step_impl( system , in , m_dxdt , t , out , m_dxdt , dt , xerr );
	}


    /*
     * version 8 : do_step( sys , in , dxdt_in , t , out , dxdt_out , dt )
     *
     * this version does not solve the forwarding problem, boost.range can not be used
     */
	template< class System , class StateIn , class DerivIn , class StateOut , class DerivOut , class Err >
	void do_step( System system , const StateIn &in , const DerivIn &dxdt_in , const time_type &t ,
			StateOut &out , DerivOut &dxdt_out , const time_type &dt , Err &xerr )
	{
		m_first_call = true;
		this->stepper().do_step_impl( system , in , dxdt_in , t , out , dxdt_out , dt , xerr );
	}





	template< class StateType >
	void adjust_size( const StateType &x )
	{
		if( m_size_adjuster.adjust_size( x ) )
		    m_first_call = true;
	}


private:

	template< class System , class StateInOut >
	void do_step_v1( System system , StateInOut &x , const time_type &t , const time_type &dt )
	{
		if( m_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() ) || m_first_call )
	    {
			typename boost::unwrap_reference< System >::type &sys = system;
	        sys( x , m_dxdt ,t );
	        m_first_call = false;
	    }
		this->stepper().do_step_impl( system , x , m_dxdt , t , x , m_dxdt , dt );
	}

	template< class System , class StateInOut , class Err >
	void do_step_v5( System system , StateInOut &x , const time_type &t , const time_type &dt , Err &xerr )
	{
	    if( m_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() ) || m_first_call )
	    {
	    	typename boost::unwrap_reference< System >::type &sys = system;
	        sys( x , m_dxdt ,t );
	        m_first_call = false;
	    }
		this->stepper().do_step_impl( system , x , m_dxdt , t , x , m_dxdt , dt , xerr );
	}



    stepper_type& stepper( void )
    {
    	return *static_cast< stepper_type* >( this );
    }

    const stepper_type& stepper( void ) const
    {
    	return *static_cast< const stepper_type* >( this );
    }


	size_adjuster< deriv_type , 1 > m_size_adjuster;
	deriv_type m_dxdt;
	bool m_first_call;

};


} // odeint
} // numeric
} // boost

#endif //BOOST_NUMERIC_ODEINT_EXPLICIT_STEPPER_AND_ERROR_STEPPER_FSAL_BASE_HPP_INCLUDED
