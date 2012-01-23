/*
 [auto_generated]
 boost/numeric/odeint/stepper/controller/pi_controller.hpp

 [begin_description]
 PI-controller for the use in the generic controller steppers.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLER_PI_CONTROLLER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLER_PI_CONTROLLER_HPP_INCLUDED

#include <cmath>

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>

namespace boost {
namespace numeric {
namespace odeint {


template< class Value >
class pi_controller
{
public:

    typedef Value value_type;

    pi_controller( void )
    : m_err_old( static_cast< value_type >( 1.0e-4 ) ) , m_last_rejected( false ) { }


    /*
     * Returns true if err Ä 1, false otherwise. If step was successful, sets hnext to the estimated
     * optimal stepsize for the next step. If the step failed, reduces h appropriately for another try.
     *
     * Set beta to a nonzero value for PI control. beta D 0:04–0.08 is a good default.
     */
    template< class Time , class Order >
    controlled_step_result operator()( value_type error , Time &t , Time &dt , Order stepper_order , Order error_order )
    {
        using std::max;
        using std::min;
        using std::pow;

        // ToDo : adapt for values of stepper_order or error_order
        const value_type beta = 0.02;
        const value_type alpha = 0.2 - beta * 0.75;
        const value_type safe = 0.9;
        const value_type min_scale = 0.2;
        const value_type max_scale = 10.0;



        if( error <= 1.0 )
        {
            value_type scale = 1.0;
            if( error == 0.0 )
            {
                scale = max_scale;
            }
            else
            {
                // PI control if beta != 0.
                scale = safe * pow( error ,-alpha ) * pow( m_err_old , beta );
                if( scale < min_scale )
                    scale = min_scale;
                if( scale > max_scale )
                    scale = max_scale;
            }

            if( m_last_rejected )
            {
                dt *= min( scale , 1.0 );
            }
            else
            {
                dt *= scale;
            }
            m_err_old = std::max( error , 1.0e-4 );
            m_last_rejected = false;
            t += dt;
            return success;
        }
        else
        {
            value_type scale = max( safe * pow( error , -alpha ) , min_scale );
            dt *= scale;
            m_last_rejected = true;
            return fail;
        }
    }

private:

    value_type m_err_old;
    bool m_last_rejected;
};



}
}
}

#endif // BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLER_PI_CONTROLLER_HPP_INCLUDED
