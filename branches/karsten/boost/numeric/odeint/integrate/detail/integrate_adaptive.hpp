/*
 * integrate_adaptive.hpp
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_ADAPTIVE_HPP_
#define BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_ADAPTIVE_HPP_

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {



template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_adaptive(
		Stepper stepper , System system , State &start_state ,
		Time &start_time , const Time &end_time , Time &dt ,
		Observer &observer , stepper_tag
		)
{
	size_t count = 0;
	while( start_time < end_time )
	{
		stepper.do_step( system , start_state , start_time , dt );
		observer( start_time , start_state );
		start_time += dt;
		++count;
	}
	return count;
}




template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_adaptive(
		Stepper stepper , System system , State &start_state ,
		Time &start_time , const Time &end_time , Time &dt ,
		Observer observer , controlled_stepper_tag
		)
{
	size_t num_of_steps = 0;
	while( start_time < end_time )
	{
		observer( start_time , start_state );
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
		while( ( res == step_size_decreased ) && ( trials < 1000 ) );
		++num_of_steps;
	}
	observer( start_time , start_state );
	return num_of_steps;
}



template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_adaptive(
		Stepper stepper , System system , State &start_state ,
		const Time &start_time , const Time &end_time , const Time &dt ,
		Observer observer , dense_output_stepper_tag )
{
	return 0;
}




} // namespace detail
} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_ADAPTIVE_HPP_ */
