/*
 * nr_rk4.cpp
 *
 *  Created on: Apr 28, 2011
 *      Author: mario
 */

#include <iostream>
#include <fstream>

#include <boost/array.hpp>

#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>
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

struct lorenz {

    static const double sigma = 10.0;
    static const double R = 28.0;
    static const double b = 8.0 / 3.0;

    void operator()( const state_type &x , state_type &dxdt , const double t ) const
    {
        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = x[0]*x[1] - b * x[2];
    }
};


template< class System , typename T , size_t dim >
void rk4_step( System sys , boost::array< T , dim > &x , const double t , const double dt )
{   // fast rk4 implementation adapted from the book 'Numerical Recipes'
    size_t i;
    const double hh = dt*0.5;
    const double h6 = dt/6.0;
    const double th = t+hh;
    boost::array< T , dim > dydx , dym , dyt , yt;

    sys( x , dydx , t );

    for( i=0 ; i<dim ; i++ )
        yt[i] = x[i] + hh*dydx[i];

    sys( yt , dyt , th );
    for( i=0 ; i<dim ; i++ )
        yt[i] = x[i] + hh*dyt[i];

    sys( yt , dym , th );
    for( i=0 ; i<dim ; i++ ) {
        yt[i] = x[i] + dt*dym[i];
        dym[i] += dyt[i];
    }
    sys( yt , dyt , t+dt );
    for( i=0 ; i<dim ; i++ )
        x[i] += h6*( dydx[i] + dyt[i] + 2.0*dym[i] );
}


const size_t loops = 20;

int main( int argc , char **argv )
{
    const size_t num_of_steps = 20000000;
    const double dt = 1E-10;

    accumulator_type acc;
    timer_type timer;

    srand( 12312354 );

    for( size_t n=0 ; n<loops ; ++n )
    {
        state_type x = {{ 10.0 * rand()/RAND_MAX ,
                          10.0 * rand()/RAND_MAX ,
                          10.0 * rand()/RAND_MAX }};
        double t = 0.0;

        timer.restart();
        for( size_t i=0 ; i<num_of_steps ; ++i, t+=dt )
            rk4_step( lorenz() , x , t , dt );
        acc( timer.elapsed() );

        clog.precision( 3 );
        clog.width( 5 );
        clog << acc << " " << x[0] << endl;
    }
    cout << acc << endl;
    return 0;
}
