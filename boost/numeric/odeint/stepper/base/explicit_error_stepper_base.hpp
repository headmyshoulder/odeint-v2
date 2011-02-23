/*
 boost header: numeric/odeint/explicit_error_stepper_base.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_EXPLICIT_ERROR_STEPPER_BASE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_EXPLICIT_ERROR_STEPPER_BASE_HPP_INCLUDED

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
 * base class for explicit error steppers
 * models the error stepper concept
 * ToDo : test
 */
template<
	class ErrorStepper ,
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
class explicit_error_stepper_base
{
public:

	typedef State state_type;
	typedef Value value_type;
	typedef Deriv deriv_type;
	typedef Time time_type;
	typedef Algebra algebra_type;
	typedef Operations operations_type;
	typedef AdjustSizePolicy adjust_size_policy;
	typedef ErrorStepper stepper_type;
	typedef explicit_error_stepper_tag stepper_category;

	typedef unsigned short order_type;
	static const order_type stepper_order_value = StepperOrder;
	static const order_type error_order_value = ErrorOrder;

    order_type stepper_order( void ) const
    {
    	return stepper_order_value;
    }

    order_type error_order( void ) const
    {
    	return error_order_value;
    }




	explicit_error_stepper_base( void ) : m_size_adjuster() , m_dxdt()
	{
		boost::numeric::odeint::construct( m_dxdt );
		m_size_adjuster.register_state( 0 , m_dxdt );
	}

	~explicit_error_stepper_base( void )
	{
		boost::numeric::odeint::destruct( m_dxdt );
	}

	explicit_error_stepper_base( const explicit_error_stepper_base &b ) : m_size_adjuster() , m_dxdt()
	{
		boost::numeric::odeint::construct( m_dxdt );
		m_size_adjuster.register_state( 0 , m_dxdt );
		boost::numeric::odeint::copy( b.m_dxdt , m_dxdt );
	}

	explicit_error_stepper_base& operator=( const explicit_error_stepper_base &b )
	{
		boost::numeric::odeint::copy( b.m_dxdt , m_dxdt );
		return *this;
	}



	/*
	 * Version 1 : do_step( system , x , t , dt , xerr )
	 *
	 * Solves for the forwading problem
	 */
	template< class System , class StateInOut , class Err >
	void do_step( System system , StateInOut &x , const time_type &t , const time_type &dt , Err &xerr )
	{
		do_step_v1( system , x , t , dt , xerr );
	}

	template< class System , class StateInOut , class Err >
	void do_step( System system , StateInOut &x , const time_type &t , const time_type &dt , Err &xerr )
	{
		do_step_v1( system , x , t , dt , xerr );
	}
	/*
	 * Version 2 : do_step( system , x , dxdt , t , dt , xerr )
	 *
	 * Does not solve for the forwarding problem
	 */
	template< class System , class StateInOut , class DerivIn , class Err >
	void do_step( System system , StateInOut &x , const DerivIn &dxdt , const time_type &t , const time_type &dt , Err &xerr )
	{
		this->stepper().do_step_impl( system , x , dxdt , t , x , dt , xerr );
	}

	/*
	 * Version 2 : do_step( system , in , t , out , dt , xerr )
	 *
	 * Does not solve for the forwarding problem
	 */
	template< class System , class StateIn , class StateOut , class Err >
	void do_step( System system , const StateIn &in , const time_type &t , StateOut &out , const time_type &dt , Err &xerr )
	{
		typename boost::unwrap_reference< System >::type &sys = system;
		m_size_adjuster.adjust_size_by_policy( in , adjust_size_policy() );
		sys( in , m_dxdt ,t );
		this->stepper().do_step_impl( system , in , m_dxdt , t , out , dt , xerr );
	}

	/*
	 * Version 2 : do_step( system , in , dxdt , t , out , dt , xerr )
	 *
	 * Does not solve for the forwarding problem
	 */
	template< class System , class StateIn , class DerivIn , class StateOut , class Err >
	void do_step( System system , const StateIn &in , const DerivIn &dxdt , const time_type &t , StateOut &out , const time_type &dt , Err &xerr )
	{
		this->stepper().do_step_impl( system , in , dxdt , t , out , dt , xerr );
	}



	template< class StateType >
	void adjust_size( const StateType &x )
	{
		m_size_adjuster.adjust_size( x );
	}


private:


	template< class System , class StateInOut , class Err >
	void do_step_v1( System system , StateInOut &x , const time_type &t , const time_type &dt , Err &xerr )
	{
		typename boost::unwrap_reference< System >::type &sys = system;
		m_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
		sys( x , m_dxdt ,t );
		this->stepper().do_step_impl( system , x , m_dxdt , t , x , dt , xerr );
	}


	// ToDo : make the next two methods private?
    stepper_type& stepper( void )
    {
    	return *static_cast< stepper_type* >( this );
    }

    const stepper_type& stepper( void ) const
    {
    	return *static_cast< const stepper_type* >( this );
    }

    void


	size_adjuster< deriv_type , 1 > m_size_adjuster;
	deriv_type m_dxdt;
};


} // odeint
} // numeric
} // boost

#endif //BOOST_NUMERIC_ODEINT_EXPLICIT_ERROR_STEPPER_BASE_HPP_INCLUDED
