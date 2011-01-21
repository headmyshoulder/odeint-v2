/*
 * check_dense_output_dopri5.cpp
 *
 *  Created on: Nov 1, 2010
 *      Author: karsten
 */

#define BOOST_TEST_MODULE odeint_dense_output_dopri5

#include <tr1/array>
#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/implicit_euler.hpp>

using namespace boost::unit_test;
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


BOOST_AUTO_TEST_SUITE( dense_output_dopri5_test )

BOOST_AUTO_TEST_CASE( test_dopri5 )
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

//    // compare with analytic solution of above system
//    BOOST_CHECK_MESSAGE( abs( x(0) - 20.0/81.0 ) < eps , x(0) - 20.0/81.0 );
//    BOOST_CHECK_MESSAGE( abs( x(1) - 10.0/9.0 ) < eps , x(0) - 10.0/9.0 );

}

BOOST_AUTO_TEST_SUITE_END()
