/*
 * test_ranges.cpp
 *
 *  Created on: Jan 18, 2011
 *      Author: karsten
 */

#include <iostream>
#include <vector>
#include <tr1/array>

#include "explicit_euler_units.hpp"

typedef std::vector< double > deriv_type;

typedef boost::numeric::odeint::explicit_euler_units< deriv_type > stepper_type;

struct harm_osc
{
	template< class State , class Deriv >
	void operator()( const State &x , Deriv &dxdt , double t )
	{
		dxdt[0] = x[1];
		dxdt[1] = -x[0];
	}
};

int main( int argc , char **argv )
{
	stepper_type stepper;
	std::tr1::array< double , 2 > x = {{ 1.0 , 0.0 }};
	stepper.do_step( harm_osc() , x , 0.0 , 0.1 );
}
