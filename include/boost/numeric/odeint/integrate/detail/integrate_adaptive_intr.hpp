/*
 [auto_generated]
 boost/numeric/odeint/integrate/detail/integrate_adaptive.hpp

 [begin_description]
 Default Integrate adaptive implementation.
 [end_description]

 Copyright 2011-2013 Karsten Ahnert
 Copyright 2011-2012 Mario Mulansky
 Copyright 2012 Christoph Koke

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_ADAPTIVE_INTR_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_ADAPTIVE_INTR_HPP_INCLUDED

#include <stdexcept>

#include <boost/throw_exception.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_const.hpp>
#include <boost/numeric/odeint/util/bind.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/numeric/odeint/util/copy.hpp>

#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>


#include <iostream>

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< class Stepper , class System , class State , class Time , class Interrupt , class Observer>
size_t integrate_adaptive_interruptible(
        Stepper stepper , System system , State &start_state ,
        Time &start_time , Time end_time , Time &dt , Interrupt interrupt,
        Observer observer , controlled_stepper_tag
)
{
    typename odeint::unwrap_reference< Observer >::type &obs = observer;
    typename odeint::unwrap_reference< Interrupt >::type &intr = interrupt;
    typename odeint::unwrap_reference< Stepper >::type &st = stepper;

    const size_t max_attempts = 1000;
    const char *error_string = "Integrate adaptive : Maximal number of iterations reached. A step size could not be found.";
    size_t count = 0;
    Time next_interrupt_time=intr(start_state,start_time);
    while( less_with_sign( start_time , end_time , dt ) )
    {
        if(start_time==next_interrupt_time)
            next_interrupt_time=intr(start_state,start_time);
        obs( start_state , start_time );
        Time force_point=(end_time<next_interrupt_time?end_time:next_interrupt_time); // minimum
        if( less_with_sign( force_point , static_cast<Time>(start_time + dt) , dt ) )
        {
            dt = force_point - start_time;
        }

        size_t trials = 0;
        controlled_step_result res = success;
        do
        {
            res = st.try_step( system , start_state , start_time , dt );
            ++trials;
        }
        while( ( res == fail ) && ( trials < max_attempts ) );
        if( trials == max_attempts ) BOOST_THROW_EXCEPTION( std::overflow_error( error_string ) );

        ++count;
    }
    obs( start_state , start_time );
    return count;
}


} // namespace detail
} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_ADAPTIVE_HPP_INCLUDED
