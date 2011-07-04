/*
 * dense_output_explicit.cpp
 *
 *  Created on: Jan 28, 2011
 *      Author: karsten
 */

#include <iostream>
#include <fstream>
#include <vector>

#include <boost/numeric/odeint/stepper/explicit_euler.hpp>
#include <boost/numeric/odeint/stepper/dense_output_explicit.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef std::vector< double > state_type;
typedef explicit_euler< state_type > explicit_stepper_type;
typedef dense_output_explicit< explicit_stepper_type > stepper_type;

struct lorenz
{
	template< class State , class Deriv >
	void operator()( const State &x , Deriv &dxdt , double t )
	{
		const double sigma = 10.0;
		const double R = 28.0;
		const double b = 8.0 / 3.0;

		dxdt[0] = sigma * ( x[1] - x[0] );
		dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
		dxdt[2] = x[0]*x[1] - b * x[2];
	}
};

int main( int argc , char **argv )
{
	stepper_type stepper;

	state_type x_start( 3 );
	x_start[0] = 10.0 , x_start[1] = 10.0 ; x_start[2] = 20.0;

	const double dt1 = 0.025 , dt2 = 0.01;
	stepper.initialize( x_start , 0.0 , dt1 );

	ofstream fout( "test.dat" );
	ofstream fout2( "test2.dat" );
	state_type x( 3 );
	double t = 0.0;
	while( stepper.current_time() < 10.0 )
	{
		while( t < stepper.current_time() )
		{
			stepper.calc_state( t , x );
			fout << t << " " << x[0] << " " << x[1] << " " << x[2] << endl;
			t += dt2;
		}
		stepper.do_step( lorenz() );
		const state_type &current = stepper.current_state();
		fout2 << stepper.current_time() << " " << current[0] << " " << current[1] << " " << current[2] << " " << endl;
	}


	return 0;
}
