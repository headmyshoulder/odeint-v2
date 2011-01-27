/* Boost check_implicit_euler.cpp test file

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 This file tests the use of the euler stepper

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#include <tr1/array>
#include <fstream>
#include <iostream>

#include <boost/numeric/odeint/stepper/explicit_euler.hpp>
#include <boost/numeric/odeint/stepper/dense_output_explicit_euler.hpp>


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

    dense_output_explicit_euler< state_type > stepper;
    state_type x0;
    x0[0] = 0.0;
    x0[1] = 1.0;

    stepper.initialize( x0 , 0.0 , 0.1 );
//    stepper.do_step( sys );

    std::ofstream stepper_out( "stepper_states.dat" );
    std::ofstream states_out( "states.dat" );


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

