/*
 * dense_output_explicit_controlled_fsal.cpp
 *
 *  Created on: Jan 28, 2011
 *      Author: karsten
 */

#include <iostream>

#include <boost/numeric/odeint/stepper/explicit_error_dopri5.hpp>
#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>
#include <boost/numeric/odeint/stepper/dense_output_controlled_explicit_fsal.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef vector< double > state_type;
typedef explicit_error_dopri5< state_type > dopri5_type;
typedef controlled_error_stepper< dopri5_type > controlled_error_stepper_type;
typedef dense_output_controlled_explicit_fsal< controlled_error_stepper_type > stepper_type;


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
	return 0;
}
