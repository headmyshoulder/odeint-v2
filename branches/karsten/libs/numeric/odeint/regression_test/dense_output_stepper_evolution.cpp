/*
 * dense_output_stepper_evolution.cpp
 *
 *  Created on: Nov 2, 2010
 *      Author: karsten
 */

#include <string>
#include <fstream>
#include <iostream>
#include <tr1/array>

#include <boost/numeric/odeint/stepper/dense_output_dopri5.hpp>
#include <boost/numeric/odeint/stepper/dense_output_explicit_euler.hpp>
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
void evolution( Stepper &stepper , double t_end , const state_type &x0 , const string &stepper_filename , const string &state_filename )
{
	ofstream state_out( state_filename.c_str() );
	ofstream stepper_out( stepper_filename.c_str() );

	stepper.initialize( x0 , 0.0 , 0.01 );

	state_type x = x0;
	double t = 0.0;
	double dt = 0.001;
	while( t < t_end )
	{
    	if( t < stepper.current_time() )
    	{
    		stepper.calc_state( t , x );
    		state_type orig = sys_solution( t , x0 );
    		state_type diff = {{ orig[0] - x[0] , orig[1] - x[1] }};
    		double diff_abs = sqrt( diff[0] * diff[0] + diff[1] * diff[1] );
    		state_out << t << tab << x << tab << orig << tab << diff << tab << diff_abs << endl;
    	}
    	else
    	{
    		stepper.do_step( sys );
    		state_type orig = sys_solution( stepper.current_time() , x0 );
    		const state_type &xx = stepper.current_state();
    		state_type diff = {{ orig[0] - xx[0] , orig[1] - xx[1] }};
    		double diff_abs = sqrt( diff[0] * diff[0] + diff[1] * diff[1] );
    		stepper_out << stepper.current_time() << "\t" << xx << tab << orig << tab << diff << tab << diff_abs << std::endl;
    		continue;
    	}
    	t += dt;

	}
}


int main( int argc , char **argv )
{
	typedef explicit_error_dopri5< state_type > dopri5_type;
	typedef controlled_error_stepper< dopri5_type > controlled_dopri5_type;

	dopri5_type dopri5;
	controlled_dopri5_type controlled_dopri5( dopri5 );

	dense_output_explicit_euler< state_type > dense_euler;
	dense_output_dopri5< controlled_dopri5_type > dense_dopri5( controlled_dopri5 );

	state_type x0 = {{ 1.25 , 0.43 }};
	evolution( dense_euler , 100.0 , x0 , "dense_euler_stepper.dat" , "dense_euler_state.dat" );
	evolution( dense_dopri5 , 100.0 , x0 , "dense_dopri5_stepper.dat" , "dense_dopri5_state.dat" );

}
