/*
 * integrate_const.hpp
 *
 *
 * Overview:
 *
 * size_t integrate( stepper , system , start_state , start_time , end_time , dt , observer );
 *
 * Time integrate_n_steps( stepper , system , start_state , start_time , dt , num_of_steps , observer );
 *
 * size_t integrate_adaptive( stepper , system , start_state , start_time , end_time , start_dt , observer );
 *
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_HPP_
#define BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_HPP_

#include <boost/numeric/odeint/stepper/explicit_error_rk54_ck.hpp>
#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>
#include <boost/numeric/odeint/integrate/do_nothing_observer.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>



namespace boost {
namespace numeric {
namespace odeint {


/*
 * ToDo :
 *
 * determine type of dxdt for units
 */
template< class System , class State , class Time , class Observer >
size_t integrate( System system , State &start_state , Time start_time , Time end_time , Time dt , Observer observer )
{
	return integrate_adaptive( controlled_error_stepper< explicit_error_rk54_ck< State > >() , system , start_state , start_time , end_time , dt , observer );
}

template< class System , class State , class Time >
size_t integrate( System system , State &start_state , Time start_time , Time end_time , Time dt )
{
	return integrate( system , start_state , start_time , end_time , dt , do_nothing_observer() );
}



} // namespace odeint
} // namespace numeric
} // namespace boost



#endif /* BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_HPP_ */
