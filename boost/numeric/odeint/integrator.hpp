/* Boost odeint/integrator.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file includes standard integration methods with adaptive stepsize.

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATOR_HPP
#define BOOST_NUMERIC_ODEINT_INTEGRATOR_HPP

#include <boost/numeric/odeint/stepsize_controller_standard.hpp>
#include <boost/numeric/odeint/resizer.hpp>
#include <vector>

namespace boost {
namespace numeric {
namespace odeint {

    class integrator {

    public:
	template< class StepType,
		  class DynamicalSystem,
		  class StateType,
		  class T >
	size_t integrate(StepType &stepper, DynamicalSystem &system, StateType &x, 
			 std::vector<T> &times, std::vector<StateType> &x_vec,
			 T dt = 1E-4, T eps_abs = 1E-7, 
			 T eps_rel = 1E-8, T a_x = 1.0 , T a_dxdt = 1.0)
	{
	    if( times.size() != x_vec.size() ) throw;
	    step_controller_standard< StateType, T >
		controller(eps_abs, eps_rel, a_x, a_dxdt );

	    typename std::vector<T>::iterator t_iter = times.begin();
	    typename std::vector<StateType>::iterator x_iter = x_vec.begin();
	    controlled_step_result result;
	    T t = *t_iter;

	    while( t_iter < times.end() ) {

		if( t >= *t_iter ) {
		    *x_iter++ = x;
		    t_iter++;
		}

		result = controller.controlled_step( stepper, system, x, t, dt );
		while( result != SUCCESS ) {
		    result = controller.controlled_step( stepper, system, x, t, dt );
		    if( dt < 1E-10 ) throw;
		}
	    }
	    return 0;
	}
    };


} // namespace odeint
} // namespace numeric
} // namespace boost

#endif
