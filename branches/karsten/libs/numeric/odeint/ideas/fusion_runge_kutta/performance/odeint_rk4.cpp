/*
 * odeint_rk4.cpp
 *
 *  Created on: Apr 28, 2011
 *      Author: mario
 */

#include <iostream>
#include <fstream>

#include <boost/array.hpp>

#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>
#include <boost/numeric/odeint/algebra/array_algebra.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/timer.hpp>

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
typedef boost::numeric::odeint::explicit_rk4< state_type > rk4_odeint_type;
//typedef boost::numeric::odeint::explicit_rk4< state_type , double , state_type , double ,
//                                              boost::numeric::odeint::array_algebra > rk4_odeint_type;


void lorenz( const state_type &x , state_type &dxdt , double t )
{
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}




int main( int argc , char **argv )
{
    rk4_odeint_type rk4_odeint;

    const size_t num_of_steps = 20000000;
    const double dt = 1E-10;

    accumulator_type acc;
    timer_type timer;

    srand( 12312354 );

    while( true )
    {
        state_type x = {{ 10.0 * rand()/RAND_MAX , 
                          10.0 * rand()/RAND_MAX , 
                          10.0 * rand()/RAND_MAX }};
        double t = 0.0;

        timer.restart();
        for( size_t i=0 ; i<num_of_steps ; ++i, t+=dt )
            rk4_odeint.do_step( lorenz , x , t , dt );
        acc( timer.elapsed() );

        clog.precision( 3 );
        clog.width( 5 );
        clog << acc << " " << x[0] << endl;
    }

    return 0;
}
