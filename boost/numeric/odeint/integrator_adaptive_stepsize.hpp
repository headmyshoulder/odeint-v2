/* Boost odeint/integrator_adaptive.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file includes integration methods with adaptive stepsize.

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATOR_ADAPTIVE_STEPSIZE_HPP
#define BOOST_NUMERIC_ODEINT_INTEGRATOR_ADAPTIVE_STEPSIZE_HPP

#include <boost/numeric/odeint/stepsize_controller_standard.hpp>
#include <boost/numeric/odeint/resizer.hpp>
#include <boost/numeric/odeint/observer.hpp>
#include <vector>
#include <limits>

namespace boost {
namespace numeric {
namespace odeint {


    template<
            class Stepper,
            class DynamicalSystem,
            class StepController,
            class Observer
            >
    size_t integrate_adaptive(
            Stepper &stepper,
            DynamicalSystem &system,
            StepController &controller,
            typename Stepper::container_type &state,
            typename Stepper::time_type start_time,
            typename Stepper::time_type end_time,
            typename Stepper::time_type dt,
            Observer &observer )
    {
        controlled_step_result result;
        size_t iterations = 0;
        typename Stepper::time_type t = start_time;
        typename Stepper::time_type dt_ = dt;

        observer(t, state, system);
        
        while( t < end_time )
        {
            // do a controlled step
            result = controller.controlled_step( stepper, system, state, t, dt_ );

            if( result != STEP_SIZE_DECREASED )
            { // we actually did a step forward (dt was small enough)
                observer(t, state, system);
                iterations++;
            }

            if( !( t+dt_ > t) ) 
                throw; // we've reached machine precision with dt - no advancing in t
        }
        return iterations;
    }

    template<
        class Stepper,
        class DynamicalSystem,
        class StepController
        >
    size_t integrate_adaptive(
            Stepper &stepper,
            DynamicalSystem &system,
            StepController &controller,
            typename Stepper::container_type &state,
            typename Stepper::time_type start_time,
            typename Stepper::time_type end_time, 
            typename Stepper::time_type dt )
    {
        return integrate_adaptive(
	    stepper , system , controller ,
	    state, start_time , end_time,
	    do_nothing_observer<
		typename Stepper::time_type ,
		typename Stepper::container_type ,
		DynamicalSystem >
	    );
    }




    /* Integration of an ode with adaptive stepsize methods.
       Integrates an ode give by system using the integration scheme stepper and the
       step-size controller controller.
       The initial state is given in x.
       t is an vector including the times at which the state will be written into 
       the vector x_vec.
       x_vec must provide enough space to hold times.size() states.
       dt is the initial step size (will be adjusted according to the controller).
       This function returns the total number of steps required to integrate the
       whole intervale times.begin() - times.end().
       Note that the values in times don't influence the stepsize, but only the 
       time points at which the state is stored into x_vec.
    */
    template< 
        class Stepper,
        class DynamicalSystem,
        class StepController,
        class TimeSequence,
        class InsertIterator
	>
    size_t integrate(
            Stepper &stepper,
            DynamicalSystem &system,
            StepController &controller,
            typename Stepper::container_type &state, 
            TimeSequence &times, 
            typename Stepper::time_type &dt, 
            InsertIterator state_inserter)
    {
        if( times.empty() ) return 0;
        else
        {
            state_copy_observer<InsertIterator, TimeSequence> observer(times, state_inserter);
            return integrate_adaptive(stepper, system, controller, state, 
                                      times.front() , times.back(), dt , observer);
        }
    }


    /* 
       Integration of an ode with adaptive stepsize methods.
       Integrates an ode give by system using the integration scheme stepper and the
       a standard step-size controller that ensures the error being below the values 
       given below.
    */
    template< 
            class Stepper,
            class DynamicalSystem,
            class InsertIterator ,
            class TimeSequence
            >
    size_t integrate(
            Stepper &stepper, 
            DynamicalSystem &system, 
            typename Stepper::container_type &x, 
            TimeSequence &times, 
            InsertIterator state_inserter,
            typename Stepper::time_type dt = 1E-4, 
            typename Stepper::time_type eps_abs = 1E-6, 
            typename Stepper::time_type eps_rel = 1E-7, 
            typename Stepper::time_type a_x = 1.0 , 
            typename Stepper::time_type a_dxdt = 1.0
                     )
    {
        // we use the standard controller for this adaptive integrator
        step_controller_standard< typename Stepper::container_type, 
	    typename Stepper::time_type, 
	    typename Stepper::resizer_type > controller(eps_abs, eps_rel, a_x, a_dxdt ); 
        // initialized with values from above
        
        // call the normal integrator
        return integrate(stepper, system, controller, x, times, dt, state_inserter);
    }
    

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif
