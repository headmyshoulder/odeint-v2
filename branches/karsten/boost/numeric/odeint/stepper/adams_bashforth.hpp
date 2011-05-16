/*
 * adams_bashforth.hpp
 *
 *  Created on: May 15, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_BASHFORTH_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_BASHFORTH_HPP_

#include <boost/ref.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/detail/adams_bashforth_coefficients.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/util/size_adjuster.hpp>


/*
 * # Introduce the number of states
 */

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

	}

public :

	typedef State state_type;
	typedef Value value_type;
	typedef Deriv deriv_type;
	typedef Time time_type;
	typedef Algebra algebra_type;
	typedef Operations operations_type;
	typedef AdjustSizePolicy adjust_size_policy;
	typedef explicit_stepper_tag stepper_category;

	static const size_t steps = Steps;

	typedef unsigned short order_type;
	static const order_type order_value = steps + 1;

	order_type order( void ) const { return order_value; }



	adams_bashforth( void )
	{
		initialize();
	}

	~adams_bashforth( void )
	{
	}

	adams_bashforth( const adams_bashforth &stepper )
	{
		initialize();
	}

	adams_bashforth& operator=( const adams_bashforth &stepper )
	{
		return *this;
	}

	template< class System , class StateInOut >
	void do_step( System system , StateInOut &x , const time_type &t , const time_type &dt )
	{
		// ToDo : implement
	}

	template< class System , class StateInOut >
	void do_step( System system , const StateInOut &x , const time_type &t , const time_type &dt )
	{
		// ToDo : implement
	}

	template< class System , class StateIn , class StateOut >
	void do_step( System system , const StateIn &in , const time_type &t , const StateOut &out , const time_type &dt )
	{
		// ToDo : implement
	}

	template< class System , class StateIn , class StateOut >
	void do_step( System system , const StateIn &in , const time_type &t , StateOut &out , const time_type &dt )
	{
		// ToDo : implement
	}

	template< class StateType >
	void adjust_size( const StateType &x )
	{
	}


private:

	boost::array< deriv_type , steps > m_steps_storage;
//	boost::circular_buffer< deriv_type* > m_previous_steps;
};




} // odeint
} // numeric
} // boost



#endif /* BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_BASHFORTH_HPP_ */
