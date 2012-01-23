/*
 [auto_generated]
 boost/numeric/odeint/stepper/controller/stiff_controller.hpp

 [begin_description]
 Controller for the use in the generic steppers in combination with solvers for stiff systems.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLER_STIFF_CONTROLLER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLER_STIFF_CONTROLLER_HPP_INCLUDED

#include <cmath>

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>

namespace boost {
namespace numeric {
namespace odeint {



template< class Value >
class stiff_controller
{
public:

    typedef Value value_type;

    stiff_controller( void )
    : m_first_step( true ) , m_last_rejected( false ) ,
      m_err_old( 0.0 ) , m_dt_old( 0.0 )
    { }

    template< class Order >
    controlled_step_result operator()( value_type error , value_type &t , value_type &dt , Order stepper_order , Order error_order )
    {
        typedef Value value_type;

        static const value_type safe = 0.9 , fac1 = 5.0 , fac2 = 1.0 / 6.0;

        value_type fac = std::max( fac2 ,std::min( fac1 , std::pow( error , 0.25 ) / safe ) );
        value_type dt_new = dt / fac;
        if ( error <= 1.0 )
        {
            if( m_first_step )
            {
                m_first_step = false;
            }
            else
            {
                value_type fac_pred = ( m_dt_old / dt ) * pow( error * error / m_err_old , 0.25 ) / safe;
                fac_pred = std::max( fac2 , std::min( fac1 , fac_pred ) );
                fac = std::max( fac , fac_pred );
                dt_new = dt / fac;
            }

            m_dt_old = dt;
            m_err_old = std::max( 0.01 , error );
            if( m_last_rejected )
                dt_new = ( dt >= 0.0 ? std::min( dt_new , dt ) : std::max( dt_new , dt ) );
            t += dt;
            dt = dt_new;
            m_last_rejected = false;
            return success;
        }
        else
        {
            dt = dt_new;
            m_last_rejected = true;
            return fail;
        }


    }

private:

    bool m_first_step;
    bool m_last_rejected;
    value_type m_err_old;
    value_type m_dt_old;
};



}
}
}

#endif // BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLER_STIFF_CONTROLLER_HPP_INCLUDED
