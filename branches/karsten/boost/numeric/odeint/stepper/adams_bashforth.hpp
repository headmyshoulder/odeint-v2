/*
 * adams_bashforth.hpp
 *
 *  Created on: May 15, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_BASHFORTH_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_BASHFORTH_HPP_

#include <boost/ref.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>

#include <boost/numeric/odeint/util/size_adjuster.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>

#include <boost/numeric/odeint/stepper/detail/adams_bashforth_coefficients.hpp>
#include <boost/numeric/odeint/stepper/detail/adams_bashforth_call_algebra.hpp>
#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>






namespace boost {
namespace numeric {
namespace odeint {


/*
 * Static explicit Adams-Bashforth multistep-solver without step size control and without dense output.
 *
 * # Define the number of steps
 */
template<
	size_t Steps ,
    class State ,
    class Value = double ,
    class Deriv = State ,
    class Time = Value ,
	class Algebra = range_algebra ,
	class Operations = default_operations ,
	class AdjustSizePolicy = adjust_size_initially_tag
	>
class adams_bashforth
{
private:

	void initialize( void )
	{
		for( size_t i=0 ; i<steps ; ++i )
			m_size_adjuster.register_state( i , m_step_storage[i] );
	}

	void copy( const adams_bashforth &stepper )
	{
		m_step_storage = stepper.m_step_storage;
	}

public :

	typedef State state_type;
	typedef Value value_type;
	typedef Deriv deriv_type;
	typedef Time time_type;
	typedef Algebra algebra_type;
	typedef Operations operations_type;
	typedef AdjustSizePolicy adjust_size_policy;
	typedef stepper_tag stepper_category;

	static const size_t steps = Steps;

	typedef unsigned short order_type;
	static const order_type order_value = steps;

	typedef detail::rotating_buffer< deriv_type , steps > step_storage_type;


	order_type order( void ) const { return order_value; }



	adams_bashforth( void )
	: m_step_storage() , m_size_adjuster() , m_coefficients()
	{
		initialize();
	}

	adams_bashforth( const adams_bashforth &stepper )
	: m_step_storage() , m_size_adjuster() , m_coefficients()
	{
		initialize();
		copy( stepper );
	}

	adams_bashforth& operator=( const adams_bashforth &stepper )
	{
		copy( stepper );
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
		m_step_storage.rotate();
		sys( in , m_step_storage[0] , t );
		detail::call_algebra< steps , algebra_type , operations_type >()( in , out , m_step_storage , m_coefficients , dt );
	}

	template< class System , class StateIn , class StateOut >
	void do_step( System system , const StateIn &in , const time_type &t , const StateOut &out , const time_type &dt )
	{
		typename boost::unwrap_reference< System >::type &sys = system;
		m_step_storage.rotate();
		sys( in , m_step_storage[0] , t );
		detail::call_algebra< steps , algebra_type , operations_type >()( in , out , m_step_storage , m_coefficients , dt );
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













	template< class StateType >
	void adjust_size( const StateType &x )
	{
		m_size_adjuster.adjust_size();
	}

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

		for( size_t i=0 ; i<steps-1 ; ++i )
		{
			if( i != 0 ) m_step_storage.rotate();
			sys( x , m_step_storage[0] , t );
			stepper.do_step( system , x , m_step_storage[0] , t , dt );
			t += dt;
		}
	}

	template< class System , class StateIn >
	void initialize( System system , StateIn &x , time_type &t , const time_type &dt )
	{
		explicit_rk4< state_type , value_type , deriv_type , time_type , algebra_type , operations_type , adjust_size_initially_tag > rk4;
		initialize( boost::ref( rk4 ) , system , x , t , dt );
	}


private:

	step_storage_type m_step_storage;
	size_adjuster< deriv_type , steps > m_size_adjuster;
	const detail::adams_bashforth_coefficients< value_type , steps > m_coefficients;
};




} // odeint
} // numeric
} // boost



#endif /* BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_BASHFORTH_HPP_ */
