/*
 * odeint_rk54ck.cpp
 *
 *  Created on: Apr 29, 2011
 *      Author: mario
 */

#include <iostream>
#include <fstream>

#include <boost/array.hpp>

#include <boost/numeric/odeint/stepper/explicit_error_rk54_ck.hpp>
#include <boost/numeric/odeint/algebra/array_algebra.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/timer.hpp>

#include "lorenz.hpp"

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


typedef boost::array< double , 3 > state_type;
//typedef boost::numeric::odeint::explicit_error_rk54_ck< state_type > rk54_ck_odeint_type;
typedef boost::numeric::odeint::explicit_error_rk54_ck< state_type , double , state_type , double ,
                                                  boost::numeric::odeint::array_algebra > rk54_ck_odeint_type;


const size_t loops = 20;

int main( int argc , char **argv )
{
    rk54_ck_odeint_type rk54_ck_odeint;

    const size_t num_of_steps = 20000000;
    double dt = 1E-10;

    accumulator_type acc;
    timer_type timer;

    srand( 12312354 );

    for( size_t n=0 ; n<loops ; ++n )
    {
        state_type x = {{ 10.0 * rand()/RAND_MAX ,
                          10.0 * rand()/RAND_MAX ,
                          10.0 * rand()/RAND_MAX }};
        state_type x_err;
        double t = 0.0;

        timer.restart();
        for( size_t i=0 ; i<num_of_steps ; ++i, t+=dt )
        {
            rk54_ck_odeint.do_step( lorenz() , x , t , dt , x_err );
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
