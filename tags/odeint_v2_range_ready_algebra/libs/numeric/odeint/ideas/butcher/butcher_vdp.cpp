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

typedef std::tr1::array< double , 2 > state_type;
typedef mpl_rk4_stepper< state_type > stepper_type;
typedef boost::numeric::odeint::explicit_rk4< state_type > stepper_std_type;

void vdp( const state_type &x , state_type &dxdt , double t )
{
	const double mu = 0.2;
	dxdt[0] = x[1];
	dxdt[1] = -x[0] + mu * ( 1.0 - x[0]*x[0] ) * x[1];
}





int main( int argc , char **argv )
{
	stepper_type midpoint;
	stepper_std_type midpoint_std;
	state_type x = {{ 1.0 , 1.0 }};
	state_type x_std = x;

	cout << "Tableau : " << endl;
	midpoint.print_tableau();
	cout << endl;

	double t = 0.0 , dt = 0.01;
	ofstream fout( "vdp.dat" );
	for( size_t i=0 ; i<10000 ; ++i )
	{
		fout << t << tab;
		fout << x[0] << tab << x[1] << tab ;
		fout << x_std[0] << tab << x_std[1] << "\n";

		midpoint.do_step( vdp , x , t , dt );
		midpoint_std.do_step( vdp , x_std , t , dt );
		t += dt;
	}



//	cout << double( one_half::num ) / double( one_half::den ) << endl;
//	cout << conv< one_half >::get_value() << endl;
//	cout << conv< one >::get_value() << endl;

	return 0;
}
