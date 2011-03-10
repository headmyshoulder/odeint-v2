/*
 * integrate_const_stepper.hpp
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_CONST_HPP_
#define BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_CONST_HPP_

#include <iostream>
using namespace std;


namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_const( Stepper stepper , System system , State &start_state , Time start_time , const Time &end_time , const Time &dt , Observer &observer , stepper_tag )
{
	while( start_time < end_time )
	{
		observer( start_time , start_state );
		stepper.do_step( system , start_state , start_time , dt );
		start_time += dt;
	}
	observer( start_time , start_state );
	return 0;
}

template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_const( Stepper stepper , System system , State &start_state , Time start_time , const Time &end_time , const Time &dt , Observer observer , controlled_stepper_tag )
{
	clog << "huhu" << endl;
	Time time_step = dt;
	while( start_time < end_time )
	{
		observer( start_time , start_state );
		Time next_time = start_time + dt;
		while( start_time < next_time )
		{
			if( ( start_time + time_step ) > next_time )
			{
				time_step = next_time - start_time;
			}
			size_t trials = 0;

			// the following loop can maybe be taken from another integrate functions
			controlled_step_result res = success_step_size_unchanged;
			do
			{
				stepper.try_step( system , start_state , start_time , time_step );
				++trials;
			}
			while( ( res == step_size_decreased ) || ( trials < 1000 ) );
		}
	}
	observer( start_time , start_state );

	return 0;
}

template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_const( Stepper stepper , System system , State &start_state , const Time &start_time , const Time &end_time , const Time &dt , Observer observer , dense_output_stepper_tag )
{
	observer( start_time , start_state );
	return 0;
}


} // namespace detail
} // namespace odeint
} // namespace numeric
} // namespace boost

#endif /* BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_CONST_HPP_ */
