/*
 * controlled_stepper_evolution.cpp
 *
 *  Created on: Nov 2, 2010
 *      Author: karsten
 */

#include <string>
#include <fstream>
#include <iostream>
#include <tr1/array>

#include <boost/numeric/odeint/stepper/explicit_error_rk54_ck.hpp>
#include <boost/numeric/odeint/stepper/explicit_error_dopri5.hpp>
#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>

#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;


typedef std::tr1::array< double , 2 > state_type;

ostream& operator<<( ostream &out , const state_type &x )
{
	out << x[0] << tab << x[1];
	return out;
}

void sys( const state_type &x , state_type &dxdt , double t )
{
	dxdt[0] = x[1];
	dxdt[1] = -x[0];
}

state_type sys_solution( double t , const state_type &x0 )
{
	state_type sol = {{ 0.0 , 0.0 }};
	sol[0] = x0[0] * cos( t ) + x0[1] * sin( t );
	sol[1] = -x0[0] * sin( t ) + x0[1] * cos( t );
	return sol;
}



template< class Stepper >
void evolution( Stepper &stepper , double t_end , const state_type &x0 , const string &filename )
{
	ofstream fout( filename.c_str() );
	state_type x = x0;
	double t = 0.0 , dt = 0.01;
	while( t < t_end )
	{
		state_type orig = sys_solution( t , x0 );
		state_type diff = {{ orig[0] - x[0] , orig[1] - x[1] }};
		double diff_abs = sqrt( diff[0] * diff[0] + diff[1] * diff[1] );
		fout << t << tab << orig << tab << x << tab << diff << tab << diff_abs << endl;
		stepper.try_step( sys , x , t , dt );
	}
}


int main( int argc , char **argv )
{
	typedef explicit_error_rk54_ck< state_type > rk54_type;
	typedef explicit_error_dopri5< state_type > dopri5_type;

	rk54_type rk54;
	dopri5_type dopri5;
	controlled_error_stepper< rk54_type > rk54_controlled( rk54 );
	controlled_error_stepper< dopri5_type > dopri5_controlled( dopri5 );

	state_type x0 = {{ 1.25 , 0.43 }};
	evolution( rk54_controlled , 100.0 , x0 , "rk54.dat" );
	evolution( dopri5_controlled , 100.0 , x0 , "dopri5.dat" );

	return 0;
}
