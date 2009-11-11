/* Boost odeint/integrator.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file includes integration methods with adaptive stepsize.

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATOR_HPP
#define BOOST_NUMERIC_ODEINT_INTEGRATOR_HPP

#include <boost/numeric/odeint/stepsize_controller_standard.hpp>
#include <boost/numeric/odeint/resizer.hpp>
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
    size_t integrate( Stepper &stepper,
                      DynamicalSystem &system,
                      StepController &controller,
                      typename Stepper::time_type start_time,
                      typename Stepper::time_type dt,
                      typename Stepper::container_type &state,
                      typename Stepper::time_type end_time,
                      Observer &observer )
    {
        controlled_step_result result;
        size_t iterations = 0;
        typename Stepper::time_type t = start_time;

        observer(t, state, system);
        
        while( t < end_time )
	{
            result = controller.controlled_step( stepper, system, state, t, dt );
            if( result != STEP_SIZE_DECREASED )
	    { // we actually did a step forward
                observer(t, state, system);
                iterations++;
            }
            while( result != SUCCESS )
	    {
                result = controller.controlled_step( stepper, system, state, t, dt );
                if( result != STEP_SIZE_DECREASED )
		{ // we did a step
                    observer(t, state, system);
                    iterations++;
                }
                if( !( t+dt > t) ) 
                    throw; // we've reached machine precision with dt - no advancing in t
            }
        }
        return iterations;
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
        class InsertIterator
	>
    size_t integrate(Stepper &stepper,
		     DynamicalSystem &system,
		     StepController &controller,
		     typename Stepper::container_type &state, 
		     std::vector<typename Stepper::time_type> &times, 
                     InsertIterator state_inserter,
		     typename Stepper::time_type dt)
    {
        // iterators for the time and state vectors
        typename std::vector<typename Stepper::time_type>::iterator t_iter = times.begin();

        controlled_step_result result;
        size_t iterations = 0;
        typename Stepper::time_type t = *t_iter;

        while( true ) { // loop will break from inside

            if( t >= *t_iter ) { // we've reached the next time point
                *state_inserter++ = state; // save the state
                t_iter++; // next time point
            }

            if( t_iter >= times.end() ) // reached end of integration time
                break; // stop loop

            result = controller.controlled_step( stepper, system, state, t, dt );
            if( result != STEP_SIZE_DECREASED )
                iterations++;
            while( result != SUCCESS ) {
                result = controller.controlled_step( stepper, system, state, t, dt );
                if( result != STEP_SIZE_DECREASED )
                    iterations++;
                if( !( t+dt > t) ) 
                    throw; // we've reached machine precision with dt - no advancing in t
            }
        }
        return iterations;
    }


    /* 
       Integration of an ode with adaptive stepsize methods.
       Integrates an ode give by system using the integration scheme stepper and the
       a standard step-size controller that ensures the error being below the values 
       given below.
       The initial state is given in x.
       t is an vector including the times at which the state will be written into 
       the vector x_vec.
       x_vec must provide enough space to hold times.size() states.
       dt is the initial step size (will be adjusted according to the errors).
       This function returns the total number of steps required to integrate the
       whole intervale times.begin() - times.end().
       Note that the values in times don't influence the stepsize, but only the 
       time points at which the state is stored into x_vec.
       
       The stepsize is adjust such that the following maximal relative error is 
       small enough for each step:
       R = max( x_err_n / [eps_abs + eps_rel*( a_x * |x_n| + a_dxdt * |dxdt_n| )] )
       where the max refers to the componentwise maximum the expression.
       
       if R > 1.1 the stepsize is decreased:
       dt = dt*S*R^(-1/q)
       
       if R < 0.5 the stepsize is increased:
       dt = dt*S*R^(-1/(q+1))

       q is the order of the stepper (e.g. 1 for simple euler) and S is a safety 
       factor set to S = 0.9.

       To avoid extensive chages in dt, the decrease factor is limited to 0.2 and 
       the increase factor to 5.0.
    */
    template< class Stepper,
              class DynamicalSystem,
              class InsertIterator
	      >
    size_t integrate(Stepper &stepper, 
                     DynamicalSystem &system, 
                     typename Stepper::container_type &x, 
                     std::vector<typename Stepper::time_type> &times, 
                     InsertIterator state_inserter,
                     typename Stepper::time_type dt = 1E-4, 
		     typename Stepper::time_type eps_abs = 1E-6, 
                     typename Stepper::time_type eps_rel = 1E-7, 
		     typename Stepper::time_type a_x = 1.0 , 
		     typename Stepper::time_type a_dxdt = 1.0)
    {
        // we use the standard controller for this adaptive integrator
        step_controller_standard< typename Stepper::container_type, 
	    typename Stepper::time_type, 
	    typename Stepper::resizer_type > controller(eps_abs, eps_rel, a_x, a_dxdt ); 
        // initialized with values from above
        
        // call the normal integrator
        return integrate(stepper, system, controller, x, times, state_inserter, dt);
    }
    

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif
