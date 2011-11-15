/*
 * stepper_details.cpp
 *
 * This example demonstrates some details about the steppers in odeint.
 *
 *  Created on: Nov 13, 2011
 *      Author: karsten
 */

#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

#include "gram_schmidt.hpp"

using namespace std;
using namespace boost::numeric::odeint;

const size_t N = 3;

typedef boost::array< double , N > state_type;

void sys( const state_type &x , state_type &dxdt , double t )
{
}






int main( int argc , char **argv )
{
    // FSAL stepper example
    {
        double t , dt;
        state_type in , out , dxdtin , dxdtout;
        //[ fsal_stepper_example
        runge_kutta_dopri5< state_type > rk;
        rk.do_step( sys , in , t , out , dt );
        rk.do_step( sys , in , dxdtin , t , dt );
        rk.do_step( sys , in , dxdtin , t , out , dxdtout , dt );
        //]
    }
    return 0;
}
