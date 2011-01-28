/*
 * dense_output_explicit.cpp
 *
 *  Created on: Jan 28, 2011
 *      Author: karsten
 */

#include <iostream>
#include <vector>

#include <boost/numeric/odeint/stepper/explicit_euler.hpp>
#include <boost/numeric/odeint/stepper/dense_output_explicit.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef std::vector< double > state_type;
typedef explicit_euler< state_type > explicit_stepper_type;
typedef dense_output_explicit< stepper_type > stepper_type;

int main( int argc , char **argv )
{
	stepper_type stepper;

	return 0;
}
