/*
 * check_dense_output_dopri5.cpp
 *
 *  Created on: Nov 1, 2010
 *      Author: karsten
 */

#include <tr1/array>
#include <fstream>
#include <iostream>

#include <boost/numeric/odeint/stepper/explicit_error_dopri5.hpp>
#include <boost/numeric/odeint/stepper/dense_output_dopri5.hpp>
#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>

using namespace boost::numeric::odeint;

typedef double value_type;
typedef std::tr1::array< double , 2 > state_type;

inline std::ostream& operator<<( std::ostream &out , const state_type &x )
{
	out << x[0] << "\t" << x[1];
	return out;
}

inline void sys( const state_type &x , state_type &dxdt , const value_type t )
{
    dxdt[0] = x[1];
    dxdt[1] = -x[0];
}


int main( int argc , char **argv )
{
    using std::abs;

    typedef explicit_error_dopri5< state_type > dopri5_type;
    typedef controlled_error_stepper< dopri5_type > controlled_stepper_type;
    dopri5_type dopri5;
    controlled_stepper_type controlled_stepper( dopri5 );
    dense_output_dopri5< controlled_stepper_type > stepper( controlled_stepper );

    state_type x0;
    x0[0] = 0.0;
    x0[1] = 1.0;

    stepper.initialize( x0 , 0.0 , 0.1 );
//    stepper.do_step( sys );

    std::ofstream stepper_out( "dopri5_stepper_states.dat" );
    std::ofstream states_out( "dopri5_states.dat" );


    double t = stepper.current_time();
    double t_end = 10.0;
    double dt = 0.02;
    state_type x;
    while( t < t_end )
    {
    	if( t < stepper.current_time() )
    	{
    		stepper.calc_state( t , x );
    		states_out << t << "\t" << x << std::endl;
    	}
    	else
    	{
    		stepper.do_step( sys );
    		stepper_out << stepper.current_time() << "\t" << stepper.current_state() << std::endl;
    		continue;
    	}
    	t += dt;
    }
}

