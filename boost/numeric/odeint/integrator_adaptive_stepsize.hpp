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

#include <boost/numeric/odeint/controlled_stepper_standard.hpp>
#include <boost/numeric/odeint/observer.hpp>
#include <vector>
#include <limits>

namespace boost {
namespace numeric {
namespace odeint {


    template<
            class ControlledStepper,
            class DynamicalSystem,
            class Observer
            >
    size_t integrate_adaptive(
            ControlledStepper &stepper,
            DynamicalSystem &system,
            typename ControlledStepper::container_type &state,
            typename ControlledStepper::time_type start_time,
            typename ControlledStepper::time_type end_time,
            typename ControlledStepper::time_type dt,
            Observer &observer )
    {
        controlled_step_result result;
        size_t iterations = 0;
        typename ControlledStepper::time_type t = start_time;
        typename ControlledStepper::time_type dt_ = dt;

        observer(t, state, system);
        
        while( t < end_time )
        {
            // do a controlled step
            result = stepper.try_step( system, state, t, dt_ );

            if( result != step_size_decreased )
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
        class ControlledStepper,
        class DynamicalSystem
        >
    size_t integrate_adaptive(
            ControlledStepper &stepper,
            DynamicalSystem &system,
            typename ControlledStepper::container_type &state,
            typename ControlledStepper::time_type start_time,
            typename ControlledStepper::time_type end_time, 
            typename ControlledStepper::time_type dt )
    {
        return integrate_adaptive(
                stepper, system ,
                state, start_time , end_time, dt, 
                do_nothing_observer<
		    typename ControlledStepper::time_type ,
		    typename ControlledStepper::container_type ,
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
        class ControlledStepper,
        class DynamicalSystem,
        class TimeInsertIterator,
        class StateInsertIterator
	>
    size_t integrate(
            ControlledStepper &stepper,
            DynamicalSystem &system,
            typename ControlledStepper::container_type &state, 
            typename ControlledStepper::time_type start,
            typename ControlledStepper::time_type end,
            typename ControlledStepper::time_type dt, 
            TimeInsertIterator time_inserter,
            StateInsertIterator state_inserter)
    {
        state_copy_observer<TimeInsertIterator, StateInsertIterator>
            observer(time_inserter, state_inserter);
        return integrate_adaptive(stepper, system, state, 
                                  start , end, dt , observer);
    }


    /* 
       Integration of an ode with adaptive stepsize methods.
       Integrates an ode give by system using the integration scheme stepper and the
       a standard step-size controller that ensures the error being below the values 
       given below.
    */
    template< 
            class DynamicalSystem,
            class ContainerType,
            class TimeInsertIterator,
            class StateInsertIterator,
            class T
            >
    size_t integrate(
            DynamicalSystem &system, 
            ContainerType &x,
            T start ,
            T end ,
            TimeInsertIterator time_inserter,
            StateInsertIterator state_inserter,
            T dt = 1E-4, 
            T eps_abs = 1E-6, 
            T eps_rel = 1E-7, 
            T a_x = 1.0 , 
            T a_dxdt = 1.0
                     )
    {
        typedef stepper_rk5_ck< ContainerType , T > stepper_type;
        // we use cash karp stepper as base stepper
        stepper_type stepper_cash_karp;
        // we use the standard controller for this adaptive integrator
        controlled_stepper_standard< stepper_type > 
            controlled_stepper(stepper_cash_karp, eps_abs, eps_rel, a_x, a_dxdt ); 
        // initialized with values from above
        
        // call the normal integrator
        return integrate(controlled_stepper, system, x, 
                         start, end, dt, time_inserter, state_inserter);
    }
    

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif
