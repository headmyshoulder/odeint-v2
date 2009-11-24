/* Boost odeint/controlled_stepper_standard.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file includes the standard controlled stepper that should be
 used with runge kutta routines.

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_CONTROLLED_STEPPER_STANDARD_HPP
#define BOOST_NUMERIC_ODEINT_CONTROLLED_STEPPER_STANDARD_HPP

#include <cmath> // for pow( ) and abs()
#include <complex>

#include <boost/concept_check.hpp>

#include <boost/numeric/odeint/detail/iterator_algebra.hpp>
#include <boost/numeric/odeint/concepts/state_concept.hpp>
#include <boost/numeric/odeint/resizer.hpp>

namespace boost {
namespace numeric {
namespace odeint {

    typedef enum{success, step_size_decreased, step_size_increased} controlled_step_result;

    /*
       The initial state is given in x.
       
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
        class ErrorStepper,
	class ResizeType = resizer< typename ErrorStepper::container_type > >
    class controlled_stepper_standard {

    public:

        // forward types from ErrorStepper
        typedef typename ErrorStepper::container_type container_type;
        typedef typename ErrorStepper::resizer_type resizer_type;
        typedef typename ErrorStepper::time_type time_type;
        typedef typename container_type::value_type value_type;
        typedef typename container_type::iterator iterator;

        typedef const unsigned short order_type;

        // private members
    private:
        ErrorStepper &m_stepper;

	time_type eps_abs;
	time_type eps_rel;
	time_type a_x;
	time_type a_dxdt;
	container_type dxdt;
	container_type x_tmp;
	container_type x_err;
	resizer_type resizer;


        // public functions
    public:
	
	controlled_stepper_standard( 
                ErrorStepper &stepper, 
                time_type abs_err, time_type rel_err, 
                time_type factor_x, time_type factor_dxdt )
	    : m_stepper(stepper), 
              eps_abs(abs_err), eps_rel(rel_err), a_x(factor_x), a_dxdt(factor_dxdt)
	{ }


        /* Tries a controlled step with the given stepsize dt. If dt is too large,
           x remains unchanged, an appropriate stepsize is assigned to dt and 
           step_size_decreased is returned. If dt is too small, the small step is 
           accomplished, a larger stepsize is assigned to dt and step_size_increased
           is returned. If dt is appropriate, the step is accomplished and success 
           is returned.
         */
	template< class DynamicalSystem >
	controlled_step_result try_step( 
                DynamicalSystem &system, 
                container_type &x, 
                time_type &t, 
                time_type &dt )
	{
	    resizer.adjust_size(x, x_err);

	    x_tmp = x; // copy current state
	    system( x, dxdt, t ); // compute dxdt
	    m_stepper.do_step(system, x, dxdt, t, dt, x_err ); // do step forward with error

	    iterator x_start = x_tmp.begin();
	    iterator dxdt_start = dxdt.begin();
	    iterator x_err_start = x_err.begin();
	    time_type max_rel_err = 0.0;

	    while( x_start != x_tmp.end() ) {
                // get the maximal value of x_err/D where 
                // D = eps_abs + eps_rel * (a_x*|x| + a_dxdt*|dxdt|);
                time_type err = eps_abs + eps_rel * (a_x * std::abs(*x_start++) + 
                                                     a_dxdt * dt * std::abs(*dxdt_start++));
                max_rel_err = max( std::abs(*x_err_start++)/err , max_rel_err );
	    }

	    //std::cout<<max_rel_err<<std::endl;

	    if( max_rel_err > 1.1 ) { // error too large - decrease dt
                // limit scaling factor to 0.2
                dt *= max( 0.9*pow(max_rel_err , -1.0/(m_stepper.order_error()-1.0)) , 0.2 );
                // reset state
                x = x_tmp;
                return step_size_decreased;
	    } else if( max_rel_err < 0.5 ) { //error too small - increase dt
                t += dt; // we keep the evolution -> increase time
                // limit scaling factor to 5.0
                dt *= min( 0.9*pow(max_rel_err , -1.0/m_stepper.order()), 5.0 );
                return step_size_increased;
	    } else {
                t += dt;
                return success;
	    }
	}
    };

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_STEPSIZE_CONTROLLER_STANDARD_HPP
