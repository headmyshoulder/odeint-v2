/*
 [auto_generated]
 boost/numeric/odeint/stepper/generic_controlled_stepper_explicit.hpp

 [begin_description]
 tba
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DEFAULT_CONTROLLER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_DEFAULT_CONTROLLER_HPP_INCLUDED

#include <cmath>

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>

namespace boost {
namespace numeric {
namespace odeint {

class default_controller
{
public:

    template< class Value , class Time , class Order >
    controlled_step_result operator()( Value error , Time &t , Time &dt , Order stepper_order , Order error_order )
    {
        using std::max;
        using std::min;
        using std::pow;

        if( error > 1.0 )
        {
            // error too large - decrease dt ,limit scaling factor to 0.2 and reset state
            dt *= max( 0.9 * pow( error , -1.0 / ( Value( error_order ) - 1.0 ) ) , 0.2 );
            return fail;
        }
        else
        {
            if( error < 0.5 )
            {
                //error too small - increase dt and keep the evolution and limit scaling factor to 5.0
                t += dt;
                dt *= min( 0.9 * pow( error , -1.0 / Value( stepper_order ) ) , 5.0 );
                return success;
            }
            else
            {
                t += dt;
                return success;
            }
        }

    }
};


}
}
}

#endif // BOOST_NUMERIC_ODEINT_STEPPER_DEFAULT_CONTROLLER_HPP_INCLUDED
