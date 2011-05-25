/*
 * adam_moulton.hpp
 *
 *  Created on: May 15, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_MOULTON_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_MOULTON_HPP_

#include <boost/ref.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>

#include <boost/numeric/odeint/util/size_adjuster.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>

#include <boost/numeric/odeint/stepper/detail/adams_bashforth_call_algebra.hpp>
#include <boost/numeric/odeint/stepper/detail/adams_moulton_coefficients.hpp>
#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>






namespace boost {
namespace numeric {
namespace odeint {


/*
 * Static implicit Adams-Moulton multistep-solver without step size control and without dense output.
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
class adam_moulton
{
private:

	void initialize( void )
	{
	}

	void copy( const adam_moulton &stepper )
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
	typedef stepper_tag stepper_category;

	static const size_t steps = Steps;

	typedef unsigned short order_type;
	static const order_type order_value = steps + 1;

	typedef detail::rotating_buffer< deriv_type , steps > step_storage_type;


	order_type order( void ) const { return order_value; }



	adam_moulton( void )
	: m_coefficients()
	{
		initialize();
	}

	adam_moulton( const adam_moulton &stepper )
	: m_coefficients()
	{
		initialize();
		copy( stepper );
	}

	adam_moulton& operator=( const adam_moulton &stepper )
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
	void do_step( System system , StateInOut &x , const time_type &t , const time_type &dt , const )
	{
		do_step( system , x , t , x , dt );
	}

	template< class System , class StateInOut >
	void do_step( System system , const StateInOut &x , const time_type &t , const time_type &dt )
	{
		do_step( system , x , t , x , dt );
	}















private:

	const detail::adams_moulton_coefficients< value_type , steps > m_coefficients;
};




} // odeint
} // numeric
} // boost



#endif /* BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_MOULTON_HPP_ */
