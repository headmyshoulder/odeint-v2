/*
 * gsl_rk4.cpp
 *
 *  Created on: Apr 28, 2011
 *      Author: mario
 */

#include <iostream>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/timer.hpp>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

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


int lorenz_gsl( const double t , const double y[] , double f[] , void *params)
{
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;

    f[0] = sigma * ( y[1] - y[0] );
    f[1] = R * y[0] - y[1] - y[0] * y[2];
    f[2] = y[0]*y[1] - b * y[2];
    return GSL_SUCCESS;
}


const size_t loops = 20;

int main()
{
    const size_t num_of_steps = 20000000 / 3 ; // gsl rk4 routine makes error control by
                                               // additional doing two steps with half step size
    const double dt = 1E-10 * 3 ;             // so it actually does 3 * num_of_steps steps

    accumulator_type acc;
    timer_type timer;

    srand( 12312354 );

    gsl_odeiv_step *s = gsl_odeiv_step_alloc( gsl_odeiv_step_rk4 , 3);
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
            gsl_odeiv_step_apply ( s , t , dt , x , x_err , 0 , 0 , &sys );
        acc( timer.elapsed() );

        clog.precision( 3 );
        clog.width( 5 );
        clog << acc << " " << x[0] << endl;
    }
    cout << acc << endl;
    gsl_odeiv_step_free( s );
    return 0;

}
