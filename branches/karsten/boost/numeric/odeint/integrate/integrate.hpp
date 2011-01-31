/*
 * integrate_const.hpp
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_HPP_
#define BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_CONST_HPP_

#include <boost/type_traits/is_same.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_const.hpp>

namespace boost {
namespace numeric {
namespace odeint {


/*

Overview:

size_t integrate( stepper , system , start_state , start_time , end_time , dt , observer );

Time integrate_n_steps( stepper , system , start_state , start_time , num_of_steps , observer );

void integrate_adaptive( stepper , system , start_state , start_time , end_time , start_dt , observer );
void integrate_adaptive( stepper , system , start_state , start_time , end_time , start_dt , state_back_insert_sequence , time_back_insert_sequence );

*/



template< class Stepper , class system , class State , class Time , class Observer = empty_observer >
size_t integrate( Stepper stepper , System system , State &start_state , const Time &start_time , const Time &end_time , const Time &dt , Observer observer = Observer() )
{
	// we want to get as fast as possible to the end
	if( typename boost:is_same< empty_observer , Observer >::value )
	{
		return detail::integrate_adpative( stepper , system , start_state , start_time , end_time  , dt , Observer , typename Stepper::stepper_type() );
	}
	else
	{
		return detail::integrate_const( stepper , system , start_state , start_time , end_time  , dt , Observer , typename Stepper::stepper_type() );
	}
}





/*
 * Old integrate functions
 */

// Constant integrator
//size_t integrate_const( stepper , system , state , start_time , end_time , dt , observer );
//size_t integrate_const( stepper , system , state , start_time , end_time , dt );
//time_type integrate_const_steps( stepper , system , state , start_time , dt , num_of_steps , observer );
//time_type integrate_const_steps( stepper , system , state , start_time , dt , num_of_steps );
//
//// Adaptive integrators
//size_t integrate_adaptive( stepper, system , state , start_time , end_time , dt , observer );
//size_t integrate_adaptive( stepper, system , state , start_time , end_time , dt );
//size_t integrate( stepper ,  system , state , start_time , end_end , dt , time_inserter , state_inserter );
//size_t integrate( system , state,start_time , end_time , time_inserter , state_inserter , dt = 1E-4 , eps_abs = 1E-6 , eps_rel = 1E-7 , a_x = 1.0 , a_dxdt = 1.0 );




} // namespace odeint
} // namespace numeric
} // namespace boost



#endif /* BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_CONST_HPP_ */
