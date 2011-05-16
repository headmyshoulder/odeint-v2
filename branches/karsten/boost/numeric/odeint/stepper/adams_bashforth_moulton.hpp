/*
 * adams_bashforth_moulton.hpp
 *
 *  Created on: May 15, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_BASHFORTH_MOULTON_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_BASHFORTH_MOULTON_HPP_

#include <boost/ref.hpp>

/*
 * # Introduce the number of states
 */

namespace boost {
namespace numeric {
namespace odeint {


/*
 * Static Adams-Bashforth-Moulton multistep-solver without step size control and without dense output.
 *
 * # Define the number of steps
 * # Define the explicit method
 * # Define the implicit method
 * # Bring the explicit method and the implicit method together
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
class adams_bashforth_moulton
{

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



	adams_bashforth_moulton( void )
	{
	}

	~adams_bashforth_moulton( void )
	{
	}

	adams_bashforth_moulton( const adams_bashforth_moulton &stepper )
	{
	}

	adams_bashforth_moulton& operator=( const adams_bashforth_moulton &stepper )
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


};




} // odeint
} // numeric
} // boost



#endif /* BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_BASHFORTH_MOULTON_HPP_ */
