/*
 * test.cpp
 *
 *  Created on: Nov 13, 2010
 *      Author: mario
 */

#include <iostream>
#include <vector>
#include <tr1/array>

#include "generic_stepper.hpp"

using namespace std;

typedef tr1::array< double , 3 > state_type;
typedef runge_kutta_stepper< state_type , 1 > euler_stepper;


const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

const double dt = 0.001;

int main( void )
{
	vector< vector <double> > a( 0, vector<double>( 0 ) );
	vector< double > b( 1 , 1.0 );
	vector< double > c( 1 , 0.0 );
	euler_stepper euler( a , b , c );

	state_type x = {{ 1.0 , 1.0 , 2.0 }};
	double t = 0.0;
	euler.do_step( lorenz , x , t , dt );

	cout << x[0] << " , " << x[1] << " , " << x[2] << endl;
}
