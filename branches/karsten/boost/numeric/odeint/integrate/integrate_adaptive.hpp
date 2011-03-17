/*
 *  integrate_adaptive.hpp
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_ADAPTIVE_HPP_
#define BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_ADAPTIVE_HPP_

#include <boost/type_traits/is_same.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/integrate/do_nothing_observer.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_const.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_adaptive.hpp>

namespace boost {
namespace numeric {
namespace odeint {



template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_adaptive(
		Stepper stepper , System system , State &start_state ,
		Time start_time , Time end_time , Time dt ,
		Observer observer )
{
	return detail::integrate_adaptive(
			stepper , system , start_state ,
			start_time , end_time , dt ,
			observer , typename Stepper::stepper_category() );
}

template< class Stepper , class System , class State , class Time >
size_t integrate_adaptive(
		Stepper stepper , System system , State &start_state ,
		Time start_time , Time end_time , Time dt )
{
	return integrate_adaptive( stepper , system , start_state , start_time , end_time , dt , do_nothing_observer() );
}




} // namespace odeint
} // namespace numeric
} // namespace boost



#endif /* BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_ADAPTIVE_HPP_ */
