/* Boost odeint/stepsize_controller_standard.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file includes the standard step size controller

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPSIZE_CONTROLLER_STANDARD_HPP
#define BOOST_NUMERIC_ODEINT_STEPSIZE_CONTROLLER_STANDARD_HPP

#include <cmath> // for pow( ) and abs()
#include <complex>

#include <boost/concept_check.hpp>

#include <boost/numeric/odeint/detail/iterator_algebra.hpp>
#include <boost/numeric/odeint/concepts/state_concept.hpp>
#include <boost/numeric/odeint/resizer.hpp>

namespace boost {
namespace numeric {
namespace odeint {

    typedef enum{SUCCESS, STEP_SIZE_DECREASED, STEP_SIZE_INCREASED} controlled_step_result;

    /*
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
       dt = dt*S*R^(-1/(q-1))
       
       if R < 0.5 the stepsize is increased:
       dt = dt*S*R^(-1/q))

       q is the order of the error term provided by the stepper.order_error() function
       (e.g. 2 for simple euler and 5 for rk5_ck) and S is a safety factor set to 
       S = 0.9. See e.g. Numerical Recipes for more details on this strategy.

       To avoid extensive chages in dt, the decrease factor is limited to 0.2 and 
       the increase factor to 5.0.
    */

    template< 
	class ContainerType, 
	class T,  
	class ResizeType = resizer< ContainerType > >
    class step_controller_standard {

	typedef typename ContainerType::iterator iterator;

	T eps_abs;
	T eps_rel;
	T a_x;
	T a_dxdt;
	ContainerType dxdt;
	ContainerType x_tmp;
	ContainerType x_err;
	ResizeType resizer;
	
    public:

	typedef ContainerType container_type;
	
	step_controller_standard( T abs_err, T rel_err, T factor_x, T factor_dxdt )
	    : eps_abs(abs_err), eps_rel(rel_err), a_x(factor_x), a_dxdt(factor_dxdt)
	{ }

	template< class Step, class DynamicalSystem >
	controlled_step_result controlled_step( Step &stepper, 
						DynamicalSystem &system, 
						ContainerType &x, 
						T &t, 
						T &dt )
	{
	    resizer.adjust_size(x, x_err);

	    x_tmp = x; // copy current state
	    system( x, dxdt, t ); // compute dxdt
	    stepper.next_step(system, x, dxdt, t, dt, x_err ); // do step forward with error

	    iterator x_start = x_tmp.begin();
	    iterator dxdt_start = dxdt.begin();
	    iterator x_err_start = x_err.begin();
	    T max_rel_err = 0.0;

	    while( x_start != x_tmp.end() ) {
                // get the maximal value of x_err/D where 
                // D = eps_abs + eps_rel * (a_x*|x| + a_dxdt*|dxdt|);
                T err = eps_abs + eps_rel * (a_x * std::abs(*x_start++) + 
                                             a_dxdt * dt * std::abs(*dxdt_start++));
                max_rel_err = max( std::abs(*x_err_start++)/err , max_rel_err );
	    }

	    //std::cout<<max_rel_err<<std::endl;

	    if( max_rel_err > 1.1 ) { // error too large - decrease dt
                // limit scaling factor to 0.2
                dt *= max( 0.9*pow(max_rel_err , -1.0/(stepper.order_error()-1.0)) , 0.2 );
                // reset state
                x = x_tmp;
                return STEP_SIZE_DECREASED;
	    } else if( max_rel_err < 0.5 ) { //error too small - increase dt
                t += dt; // we keep the evolution -> increase time
                // limit scaling factor to 5.0
                dt *= min( 0.9*pow(max_rel_err , -1.0/stepper.order()), 5.0 );
                return STEP_SIZE_INCREASED;
	    } else {
                t += dt;
                return SUCCESS;
	    }
	}
    };

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_STEPSIZE_CONTROLLER_STANDARD_HPP
