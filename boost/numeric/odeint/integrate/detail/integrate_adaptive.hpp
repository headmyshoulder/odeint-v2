/*
 * integrate_adaptive.hpp
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_ADAPTIVE_HPP_
#define BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_ADAPTIVE_HPP_

#include <stdexcept>

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {



/*
 * integrate_adaptive for simple stepper is a integrate_const
 */
template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_adaptive(
		Stepper stepper , System system , State &start_state ,
		Time start_time , Time end_time , Time dt ,
		Observer observer , stepper_tag
		)
{
	return integrate_const( stepper , system , start_state ,
			start_time , end_time , dt , observer , stepper_tag() );
}



/*
 * classical integrate adpative
 */
template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_adaptive(
		Stepper stepper , System system , State &start_state ,
		Time &start_time , Time end_time , Time &dt ,
		Observer observer , controlled_stepper_tag
		)
{
	const size_t max_attempts = 1000;
	const char *error_string = "Integrate adaptive : Maximal number of iterations reached. A step size could not be found.";
	size_t count = 0;
	while( start_time < end_time )
	{
		if( ( start_time + dt ) > end_time )
		{
			dt = end_time - start_time;
		}

		size_t trials = 0;
		controlled_step_result res = success_step_size_unchanged;
		do
		{
			res = stepper.try_step( system , start_state , start_time , dt );
			++trials;
		}
		while( ( res == step_size_decreased ) && ( trials < max_attempts ) );
		if( trials == max_attempts ) throw std::overflow_error( error_string );

		observer( start_state , start_time );
		++count;
	}
	return count;
}


/*
 * integrate adaptive for dense output steppers
 *
 * step size control is used if the stepper supports it
 */
template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_adaptive(
		Stepper stepper , System system , State &start_state ,
		Time start_time , Time end_time , Time dt ,
		Observer observer , dense_output_stepper_tag )
{
	size_t count = 0;
	stepper.initialize( start_state , start_time , dt );
	while( stepper.current_time() < end_time )
	{
		stepper.do_step( system );
		observer( stepper.current_state() , stepper.current_time() );
		++count;
	}
	return count;
}




} // namespace detail
} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_ADAPTIVE_HPP_ */
