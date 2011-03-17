/*
 * integrate_n_steps.hpp
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_N_STEPS_HPP_
#define BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_N_STEPS_HPP_

#include <boost/type_traits/is_same.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/integrate/do_nothing_observer.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_const.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_adaptive.hpp>

namespace boost {
namespace numeric {
namespace odeint {


/*
 * Integrates n steps
 */
template< class Stepper , class System , class State , class Time , class Observer >
Time integrate_n_steps(
		Stepper stepper , System system , State &start_state ,
		Time start_time , Time dt , size_t num_of_steps ,
		Observer observer )
{
	Time end_time = dt * num_of_steps;

	// we want to get as fast as possible to the end
	if( boost::is_same< do_nothing_observer , Observer >::type::value )
	{
		detail::integrate_adaptive(
				stepper , system , start_state ,
				start_time , end_time  , dt ,
				observer , typename Stepper::stepper_category() );
	}
	else
	{
		detail::integrate_const(
				stepper , system , start_state ,
				start_time , end_time  , dt ,
				observer , typename Stepper::stepper_category() );
	}
	return end_time;
}

template< class Stepper , class System , class State , class Time >
Time integrate_n_steps(
		Stepper stepper , System system , State &start_state ,
		Time start_time , Time dt , size_t num_of_steps )
{
	return integrate_n_steps( stepper , system , start_state , start_time , dt , num_of_steps , do_nothing_observer() );
}



} // namespace odeint
} // namespace numeric
} // namespace boost



#endif /* BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_N_STEPS_HPP_ */
