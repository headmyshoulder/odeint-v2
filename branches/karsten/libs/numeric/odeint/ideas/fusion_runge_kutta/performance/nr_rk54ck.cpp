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

inline void lorenz( const state_type &x , state_type &dxdt , const double t )
{
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}


template< class System , typename T , size_t dim >
void rk54ck_step( System &sys , boost::array< T , dim > &x , const double t , const double dt , boost::array< T , dim > &xerr )
{   // fast rk54ck implementation adapted from the book 'Numerical Recipes'
    size_t i;
    static const double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875 ,
            b21 = 0.2 , b31 = 3.0 / 40.0 , b32 = 9.0 / 40.0 , b41 = 0.3 , b42 =
                    -0.9 , b43 = 1.2 , b51 = -11.0 / 54.0 , b52 = 2.5 , b53 =
                    -70.0 / 27.0 , b54 = 35.0 / 27.0 , b61 = 1631.0 / 55296.0 ,
            b62 = 175.0 / 512.0 , b63 = 575.0 / 13824.0 , b64 = 44275.0
                    / 110592.0 , b65 = 253.0 / 4096.0 , c1 = 37.0 / 378.0 , c3 =
                    250.0 / 621.0 , c4 = 125.0 / 594.0 , c6 = 512.0 / 1771.0 ,
            dc1 = c1 - 2825.0 / 27648.0 , dc3 = c3 - 18575.0 / 48384.0 , dc4 =
                    c4 - 13525.0 / 55296.0 , dc5 = -277.00 / 14336.0 , dc6 = c6
                    - 0.25;
    const size_t n = dim;
    boost::array< T , dim > dydx , ak2 , ak3 , ak4 , ak5 , ak6 , ytemp;

    sys( x , dydx , t );
    for (i=0;i<n;i++)
        ytemp[i] = x[i] + b21 * dt * dydx[i];

    sys( ytemp , ak2 , t+a2*dt );
    for (i=0;i<n;i++)
        ytemp[i] = x[i] + dt*(b31*dydx[i]+b32*ak2[i]);

    sys( ytemp , ak3 , t+a3*dt );
    for (i=0;i<n;i++)
        ytemp[i] = x[i] + dt*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);

    sys( ytemp , ak4 , t+a4*dt );
    for (i=0;i<n;i++)
        ytemp[i] = x[i] + dt*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);

    sys( ytemp, ak5 , t+a5*dt );
    for (i=0;i<n;i++)
        ytemp[i] = x[i] + dt*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);

    sys( ytemp , ak6 , t+a6*dt );
    for (i=0;i<n;i++)
        x[i] += dt*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
    for (i=0;i<n;i++)
        xerr[i] = dt*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}


const size_t loops = 20;

int main( int argc , char **argv )
{
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
            rk54ck_step( lorenz , x , t , dt , x_err );
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
