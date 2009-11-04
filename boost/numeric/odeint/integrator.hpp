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
    template< class StepType,
	      class DynamicalSystem,
	      class StepController,
	      class T >
    size_t integrate(StepType &stepper,
		     DynamicalSystem &system,
		     StepController &controller,
		     typename StepType::container_type &x, 
		     std::vector<T> &times, 
		     std::vector<typename StepType::container_type> &x_vec,
		     T dt)
    {
	if( times.size() != x_vec.size() ) throw;

	// iterators for the time and state vectors
	typename std::vector<T>::iterator t_iter = times.begin();
	typename std::vector<typename StepType::container_type>::iterator x_iter = x_vec.begin();

	controlled_step_result result;
	size_t iterations = 0;
	T t = *t_iter;

	while( t_iter < times.end() ) {

	    if( t >= *t_iter ) { // we've reached the next time point
		*x_iter++ = x; // save the vector
		t_iter++; // next time point
	    }

	    result = controller.controlled_step( stepper, system, x, t, dt );
	    while( result != SUCCESS ) {
		result = controller.controlled_step( stepper, system, x, t, dt );
		if( result == STEP_SIZE_INCREASED )
		    iterations++;
		if( !( t+dt > t) ) 
		    throw; // we've reached machine precision with dt - no advancing in t
	    }
	    iterations++;
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
    template< class StepType,
	      class DynamicalSystem,
	      class T >
    size_t integrate(StepType &stepper, 
		     DynamicalSystem &system, 
		     typename StepType::container_type &x, 
		     std::vector<T> &times, 
		     std::vector<typename StepType::container_type> &x_vec,
		     T dt = 1E-4, T eps_abs = 1E-6, 
		     T eps_rel = 1E-7, T a_x = 1.0 , T a_dxdt = 1.0)
    {
	if( times.size() != x_vec.size() ) throw;
	// we use the standard controller for this adaptive integrator
	step_controller_standard< typename StepType::container_type, T, typename StepType::resizer_type>
	    controller(eps_abs, eps_rel, a_x, a_dxdt ); // initialized with values from above
	
	// call the normal integrator
	return integrate(stepper, system, controller, x, times, x_vec, dt);
    }
    

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif
