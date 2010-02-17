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
#include <algorithm>
#include <complex>

#include <boost/concept_check.hpp>

#include <boost/numeric/odeint/detail/iterator_algebra.hpp>
#include <boost/numeric/odeint/concepts/state_concept.hpp>
#include <boost/numeric/odeint/error_checker_standard.hpp>
#include <boost/numeric/odeint/container_traits.hpp>

namespace boost {
namespace numeric {
namespace odeint {

    typedef enum{
        success ,
        step_size_decreased ,
        step_size_increased
    } controlled_step_result;


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

    template< class ErrorStepper >
    class controlled_stepper_standard
    {

    public:

        // forward types from ErrorStepper
        typedef ErrorStepper stepper_type;
        typedef typename stepper_type::order_type order_type;
        typedef typename stepper_type::container_type container_type;
        typedef typename stepper_type::time_type time_type;
        typedef typename stepper_type::traits_type traits_type;
        typedef typename stepper_type::value_type value_type;
        typedef typename stepper_type::iterator iterator;
        typedef typename stepper_type::const_iterator const_iterator;


        // private members
    private:

        stepper_type &m_stepper;
        error_checker_standard< container_type, time_type , traits_type > m_error_checker;

	time_type m_eps_abs;
	time_type m_eps_rel;
	time_type m_a_x;
	time_type m_a_dxdt;
	container_type m_dxdt;
	container_type m_x_tmp;
	container_type m_x_err;
        container_type m_x_scale;


        // public functions
    public:
	
	controlled_stepper_standard( 
                ErrorStepper &stepper, 
                time_type abs_err, time_type rel_err, 
                time_type factor_x, time_type factor_dxdt )
	    : m_stepper(stepper), 
              m_error_checker( abs_err, rel_err, factor_x, factor_dxdt ),
              m_eps_abs(abs_err),
              m_eps_rel(rel_err),
              m_a_x(factor_x),
              m_a_dxdt(factor_dxdt)
	{ }

        order_type order() { return m_stepper.order(); }

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
                container_type &dxdt,
                time_type &t, 
                time_type &dt )
	{
            traits_type::adjust_size( x , m_x_err );
            traits_type::adjust_size( x , m_x_scale );

            m_error_checker.fill_scale( x , dxdt , dt , m_x_scale );

	    m_x_tmp = x;
	    m_stepper.do_step( system , x , dxdt, t , dt , m_x_err );

            time_type max_rel_err = m_error_checker.get_max_error_ratio(m_x_err, m_x_scale);

	    if( max_rel_err > 1.1 )
            { 
                // error too large - decrease dt
                // limit scaling factor to 0.2
                dt *= std::max( 0.9*pow(max_rel_err , -1.0/(m_stepper.order_error()-1.0)),
                                0.2 );

                // reset state
                x = m_x_tmp;
                return step_size_decreased;
	    }
            else
            {
                if( max_rel_err < 0.5 )
                {
                    //error too small - increase dt and keep the evolution
                    t += dt;
                    // limit scaling factor to 5.0
                    dt *= std::min( 0.9*pow(max_rel_err , -1.0/m_stepper.order()), 5.0 );
                    return step_size_increased;
                }
                else
                {
                    t += dt;
                    return success;
                }
            }
	}

	template< class DynamicalSystem >
	controlled_step_result try_step( 
                DynamicalSystem &system, 
                container_type &x, 
                time_type &t, 
                time_type &dt )
	{
            traits_type::adjust_size( x , m_dxdt );
            system( x , m_dxdt , t );
            return try_step( system , x , m_dxdt , t , dt );
        }

    };

    template< class ErrorStepper >
    controlled_stepper_standard< ErrorStepper > make_controlled_stepper_standard(
            ErrorStepper &stepper,
            ErrorStepper::time_type abs_err, ErrorStepper::time_type rel_err, 
            ErrorStepper::time_type factor_x, ErrorStepper::time_type factor_dxdt )
    {
        return controlled_stepper_standard< ErrorStepper >(
                stepper , abs_err , rel_err , factor_x , factor_dxdt );
    };
        
            

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_STEPSIZE_CONTROLLER_STANDARD_HPP
