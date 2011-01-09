/*
 * rosenbrock4.cpp
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

template< class StateType >
void system( const StateType &x , StateType &dxdt , time_type t )
{
	dxdt[0] = -0.013 * x[0] - 1000.0 * x[0] * x[2];
	dxdt[1] = -2500.0 * x[1] * x[2];
	dxdt[2] = -0.013 * x[0] - 1000.0 * x[0] * x[2] - 2500.0 * x[1] * x[2];
}

void jacobi( const state_type &x , matrix_type &J , time_type t , state_type &dfdt )
{
	J( 0 , 0 ) = -0.013 - 1000.0 * x[2];
	J( 0 , 1 ) = 0.0;
	J( 0 , 2 ) = -1000.0 * x[0];
	J( 1 , 0 ) = 0.0;
	J( 1 , 1 ) = -2500.0 * x[2];
	J( 1 , 2 ) = -2500.0 * x[1];
	J( 2 , 0 ) = -0.013 - 1000.0 * x[2];
	J( 2 , 1 ) = -2500.0 * x[2];
	J( 2 , 2 ) = -1000.0 * x[0] - 2500.0 * x[1];

	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
}


int main( int argc , char **argv )
{
	if( true )
	{
		state_type x( dim ) , xerr( dim );
		time_type t = 0.0 , dt = 0.00001;

		stepper_type stepper;
		stepper.do_step( system< state_type > , jacobi , x , t , dt , xerr );
		controlled_stepper_type controlled_stepper;

		x[0] = 1.0 ; x[1] = 1.0 ; x[2] = 0.0;

		ofstream fout( "dat/ross.dat" );
		size_t count = 0;
		while( t < 50.0 )
		{
			fout << t << "\t" << dt << "\t" << x[0] << "\t" << x[1] << "\t" << x[2] << std::endl;
			size_t trials = 0;
			while( trials < 100 )
			{
				if( controlled_stepper.try_step( system< state_type > , jacobi , x , t , dt ) !=  step_size_decreased )
					break;
				++trials;
			}
			if( trials == 100 )
			{
				cerr << "Error : stepper did not converge! " << endl;
				break;
			}
			++count;
		}
		clog << "Rosenbrock : " << count << endl;
	}





	if( true )
	{
		typedef std::tr1::array< time_type , 3 > state_type2;
		typedef explicit_error_rk54_ck< state_type2 > stepper_type2;
		typedef controlled_error_stepper< stepper_type2 > controlled_stepper_type2;
		stepper_type2 rk_stepper;
		controlled_stepper_type2 stepper( rk_stepper );

		state_type2 x = {{ 1.0 , 1.0 , 0.0 }};
		time_type t = 0.0 , dt = 0.00001;
		ofstream fout( "dat/rk.dat" );
		size_t count = 0;
		while( t < 50.0 )
		{
			fout << t << "\t" << dt << "\t" << x[0] << "\t" << x[1] << "\t" << x[2] << std::endl;
			size_t trials = 0;
			while( trials < 100 )
			{
				if( stepper.try_step( system< state_type2 > , x , t , dt ) !=  step_size_decreased )
					break;
				++trials;
			}
			if( trials == 100 )
			{
				cerr << "Error : stepper did not converge! " << endl;
				break;
			}
			++count;
		}
		clog << "RK 54 : " << count << endl;
	}

	return 0;
}
