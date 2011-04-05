/*
 * chaotic_system.cpp
 *
 * This example demonstrates how one can use odeint to determine the Lyapunov
 * exponents of a chaotic system namely the well known Lorenz system.
 *
 *  Created on: Apr 5, 2011
 *      Author: karsten
 */

#include <iostream>

#include <tr1/array>

#include <boost/numeric/odeint.hpp>

#include "gram_schmitt.hpp"

using namespace std;
using namespace boost::numeric::odeint;

const size_t n = 3;
const size_t num_of_lyap = 3;
const size_t N = n + n*num_of_lyap;

typedef std::tr1::array< double , N > state_type;
typedef std::tr1::array< double , num_of_lyap > lyap_type;
const double dt = 0.01;


struct lorenz
{
	template< class State , class Deriv >
	void operator()( const State &x_ , Deriv &dxdt_ , double t ) const
	{
		typename boost::range_iterator< const State >::type x = boost::begin( x_ );
		typename boost::range_iterator< Deriv >::type dxdt = boost::begin( dxdt_ );

		dxdt[0]=10.0*(x[1]-x[0]);
		dxdt[1]=28.0*x[0]-x[1]-x[0]*x[2];
		dxdt[2]=-2.666666666666*x[2]+x[0]*x[1];
	}
};

void lorenz_with_lyap( const state_type &x , state_type &dxdt , double t )
{
	lorenz()( x , dxdt , t );

	for( size_t l=0 ; l<num_of_lyap ; ++l )
    {
		const double *pert = x.begin() + 3 + l * 3;
		double *dpert = dxdt.begin() + 3 + l * 3;
        dpert[0] = -10.0 * pert[0] + 10.0 * pert[1];
        dpert[1] = ( 28.0 - x[2] ) * pert[0] - pert[1] - x[0] * pert[2];
        dpert[2] = x[1] * pert[0] + x[0] * pert[1] - 2.666666666666 * pert[2];
    }
}

struct streaming
{
	template< class State , class Time >
	void operator()( const State& x , const Time &t ) const
	{
		cout << t << "\t" << x[0] << "\t" << x[1] << "\t" << x[2] << "\n";
	}
};

int main( int argc , char **argv )
{
    state_type x;
    lyap_type lyap;
    explicit_rk4< state_type > rk4;

    fill( x.begin() , x.end() , 0.0 );
    x[0] = 10.0 ; x[1] = 10.0 ; x[2] = 5.0;

    const double dt = 0.1;

    // 10000 transients steps
    integrate_n_steps( rk4 , lorenz() , std::make_pair( x.begin() , x.begin() + n ) , 0.0 , dt , 10000 );

    fill( x.begin()+n , x.end() , 0.0 );
    for( size_t i=0 ; i<num_of_lyap ; ++i ) x[n+n*i+i] = 1.0;
    fill( lyap.begin() , lyap.end() , 0.0 );

    double t = 0.0;
    size_t count = 0;
    for( size_t i=0 ; i<x.size() ; ++i ) cout << "\t" << x[i];
    cout << endl;
    while( true )
    {
    	t = integrate_n_steps( rk4 , lorenz_with_lyap , x , t , dt , 100 );
        gram_schmitt( x , lyap , n , num_of_lyap );
        ++count;
        if( count == 1 )
        {
        	for( size_t i=0 ; i<x.size() ; ++i ) cout << "\t" << x[i];
        	cout << endl;
        }


        if( !(count % 1) )
        {
            cout << t;
            for( size_t i=0 ; i<num_of_lyap ; ++i )
                cout << "\t" << lyap[i] / t ;
            for( size_t i=0 ; i<x.size() ; ++i ) cout << "\t" << x[i];
            cout << endl;
        }

        if( count == 2 ) break;
    }

    return 0;
}
