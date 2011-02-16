/*
 * butcher_test.cpp
 *
 *  Created on: Nov 5, 2010
 *      Author: karsten
 */

#include <iostream>
#include <fstream>

#include <tr1/array>

#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>

#include "predefined_steppers.hpp"

#define tab "\t"

using namespace std;

typedef std::tr1::array< double , 3 > state_type;
typedef mpl_rk4_stepper< state_type > stepper_type;
typedef boost::numeric::odeint::explicit_rk4< state_type > stepper_std_type;


const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}




int main( int argc , char **argv )
{
	stepper_type stepper;
	stepper_std_type stepper_std;
	state_type x = {{ 1.0 , 1.0 , 2.0 }};
	state_type x_std = x;

	cout << "Tableau : " << endl;
	stepper.print_tableau();
	cout << endl;

	double t = 0.0 , dt = 0.01;
	ofstream fout( "lorenz.dat" );
	for( size_t i=0 ; i<10000 ; ++i )
	{
		fout << t << tab;
		fout << x[0] << tab << x[1] << tab << x[2] << tab;
		fout << x_std[0] << tab << x_std[1] << tab << x_std[2] << "\n";

		stepper.do_step( lorenz , x , t , dt );
		stepper_std.do_step( lorenz , x_std , t , dt );
		t += dt;
	}


	return 0;
}
