/*
 * rosenbrock4_stepper.cpp
 *
 *  Created on: Jan 9, 2011
 *      Author: karsten
 */

#include <iostream>
#include <fstream>
#include <tr1/array>

#include "rosenbrock4.hpp"
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const static size_t dim = 3;
typedef double time_type;
typedef rosenbrock4< time_type > stepper_type;
typedef stepper_type::state_type state_type;
typedef stepper_type::matrix_type matrix_type;
typedef rosenbrock4_controller< time_type > controlled_stepper_type;

//template< class StateType >
//void system( const StateType &x , StateType &dxdt , time_type t )
//{
//	dxdt[0] = -0.013 * x[0] - 1000.0 * x[0] * x[2];
//	dxdt[1] = -2500.0 * x[1] * x[2];
//	dxdt[2] = -0.013 * x[0] - 1000.0 * x[0] * x[2] - 2500.0 * x[1] * x[2];
//}
//
//void jacobi( const state_type &x , matrix_type &J , time_type t , state_type &dfdt )
//{
//	J( 0 , 0 ) = -0.013 - 1000.0 * x[2];
//	J( 0 , 1 ) = 0.0;
//	J( 0 , 2 ) = -1000.0 * x[0];
//	J( 1 , 0 ) = 0.0;
//	J( 1 , 1 ) = -2500.0 * x[2];
//	J( 1 , 2 ) = -2500.0 * x[1];
//	J( 2 , 0 ) = -0.013 - 1000.0 * x[2];
//	J( 2 , 1 ) = -2500.0 * x[2];
//	J( 2 , 2 ) = -1000.0 * x[0] - 2500.0 * x[1];
//
//	dfdt[0] = 0.0;
//	dfdt[1] = 0.0;
//	dfdt[2] = 0.0;
//}


const time_type sigma = 10.0;
const time_type R = 28.0;
const time_type b = 8.0 / 3.0;

template< class StateType >
void system( const StateType &x , StateType &dxdt , time_type t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0] * x[1] - b * x[2];
}

void jacobi( const state_type &x , matrix_type &J , time_type t , state_type &dfdt )
{
	J( 0 , 0 ) = -sigma;
	J( 0 , 1 ) = sigma;
	J( 0 , 2 ) = 0.0;
	J( 1 , 0 ) = R - x[2];
	J( 1 , 1 ) = -1.0;
	J( 1 , 2 ) = -x[0];
	J( 2 , 0 ) = x[1];
	J( 2 , 1 ) = x[0];
	J( 2 , 2 ) = -b;

	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
}



int main( int argc , char **argv )
{
	const double dt = 0.01;
	size_t steps = 1000;
	double x0 = -12.0 , y0 = -12.0 , z0 = 20.0;
	if( true )
	{
		state_type x( dim ) , xerr( dim );
		time_type t = 0.0;

		stepper_type stepper;
		x[0] = x0 ; x[1] = y0 ; x[2] = z0;

		ofstream fout( "dat/ross.dat" );
		fout.precision( 14 );
		size_t count = 0;
		while( count < steps )
		{
			fout << t << "\t";
			fout << x[0] << "\t" << x[1] << "\t" << x[2] << "\t";
			fout << xerr[0] << "\t" << xerr[1] << "\t" << xerr[2] << "\t";
			fout <<std::endl;

			stepper.do_step( system< state_type > , jacobi , x , t , dt , xerr );
			++count;
			t += dt;
		}
	}





	if( true )
	{
		typedef std::tr1::array< time_type , 3 > state_type2;
		typedef explicit_error_rk54_ck< state_type2 > stepper_type2;
		stepper_type2 rk_stepper;
		state_type2 x = {{ x0 , y0 , z0 }} , xerr = {{ 0.0 , 0.0 , 0.0 }};
		time_type t = 0.0;

		ofstream fout( "dat/rk.dat" );
		fout.precision( 14 );
		size_t count = 0;
		while( count < steps )
		{
			fout << t << "\t";
			fout << x[0] << "\t" << x[1] << "\t" << x[2] << "\t";
			fout << xerr[0] << "\t" << xerr[1] << "\t" << xerr[2] << "\t";
			fout <<std::endl;

			rk_stepper.do_step( system< state_type2 > , x , t , dt , xerr );
			++count;
			t += dt;
		}
	}

	return 0;
}
