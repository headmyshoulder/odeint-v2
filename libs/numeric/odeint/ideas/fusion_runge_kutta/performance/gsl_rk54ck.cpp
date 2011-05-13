/*
 * gsl_rk54ck.cpp
 *
 *  Created on: Apr 29, 2011
 *      Author: mario
 */

#include <iostream>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/timer.hpp>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "lorenz_gsl.hpp"

#define tab "\t"

using namespace std;
using namespace boost::accumulators;

typedef accumulator_set<
    double , stats< tag::mean , tag::variance >
    > accumulator_type;

ostream& operator<<( ostream& out , accumulator_type &acc )
{
    out << boost::accumulators::mean( acc ) << tab;
//    out << boost::accumulators::variance( acc ) << tab;
    return out;
}

typedef boost::timer timer_type;

const size_t loops = 20;

int main()
{
    const size_t num_of_steps = 20000000;
    double dt = 1E-10;

    accumulator_type acc;
    timer_type timer;

    srand( 12312354 );

    gsl_odeiv_step *s = gsl_odeiv_step_alloc( gsl_odeiv_step_rkck , 3);
    gsl_odeiv_system sys = { lorenz_gsl , 0 , 3 , 0 };

    for( size_t n=0 ; n<loops ; ++n )
    {
        double x[3] = { 10.0 * rand()/RAND_MAX ,
                         10.0 * rand()/RAND_MAX ,
                         10.0 * rand()/RAND_MAX };
        double x_err[3];

        double t = 0.0;

        timer.restart();
        for( size_t i=0 ; i<num_of_steps ; ++i, t+=dt )
        {
            gsl_odeiv_step_apply ( s , t , dt , x , x_err , 0 , 0 , &sys );
            if( i % 1000 == 0 ) // simulated stepsize control
                dt += (dt*1E-3*rand())/RAND_MAX - dt*5E-4;
        }
        acc( timer.elapsed() );

        clog.precision( 15 );
        clog.width( 20 );
        clog << acc << " " << x[0] << tab << " " << x_err[0] << endl;
    }
    cout << acc << endl;
    return 0;

}
