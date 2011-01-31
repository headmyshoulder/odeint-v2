/*
 * do_nothing_observer.hpp
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_DO_NOTHING_OBSERVER_HPP_
#define BOOST_NUMERIC_ODEINT_INTEGRATE_DO_NOTHING_OBSERVER_HPP_


namespace boost {
namespace numeric {
namespace odeint {

struct do_nothing_observer
{
	template< class State , class Time >
	void operator()( const State& x , const Time &t ) const
	{

	}
};

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif /* BOOST_NUMERIC_ODEINT_INTEGRATE_DO_NOTHING_OBSERVER_HPP_ */
