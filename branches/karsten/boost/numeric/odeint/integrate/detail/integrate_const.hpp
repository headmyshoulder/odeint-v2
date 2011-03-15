/*
 * integrate_const_stepper.hpp
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_CONST_HPP_
#define BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_CONST_HPP_


#include <boost/numeric/odeint/integrate/detail/integrate_adaptive.hpp>

#include <iostream>
using namespace std;


namespace boost {
namespace numeric {
namespace odeint {
namespace detail {



template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_const(
		Stepper stepper , System system , State &start_state ,
		Time start_time , Time end_time , Time dt ,
		Observer &observer , stepper_tag
		)
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
size_t integrate_const(
		Stepper stepper , System system , State &start_state ,
		Time &start_time , const Time &end_time , Time &dt ,
		Observer observer , controlled_stepper_tag
		)
{
	size_t count = 0;
	Time time_step = dt;
	while( start_time < end_time )
	{
		observer( start_time , start_state );
		Time next_time = start_time + time_step;
		if( next_time > end_time ) next_time = end_time;
		detail::integrate_adaptive(
				stepper , system , start_state , start_time , next_time , dt ,
				do_nothing_observer() , controlled_stepper_tag() );
		++count;
	}
	observer( start_time , start_state );
	return count;
}



template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_const(
		Stepper stepper , System system , State &start_state ,
		Time start_time , Time end_time , Time dt ,
		Observer observer , dense_output_stepper_tag )
{
	clog << "dense output" << endl;
	stepper.initialize( start_state , start_time , dt );

	size_t count = 0;
	size_t count2 = 0;
	while( start_time < end_time )
	{
		while( ( start_time < stepper.current_time() ) && ( start_time < end_time ) )
		{
			stepper.calc_state( start_time , start_state );
			observer( start_time , start_state );
			start_time += dt;
			++count;
		}

		// we have not reached the end, do another real step
		if( start_time < end_time )
		{
			stepper.do_step( system );
			++count2;
		}
	}
	clog << count2 << endl;
	return count;


	observer( start_time , start_state );
	return 0;
}


} // namespace detail
} // namespace odeint
} // namespace numeric
} // namespace boost

#endif /* BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_CONST_HPP_ */
