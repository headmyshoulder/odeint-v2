/*
 * main.cpp
 *
 *  Created on: Apr 2, 2011
 *      Author: karsten
 */

#include <iostream>

#include "taylor.hpp"

#include <boost/fusion/include/make_vector.hpp>

template< typename T , size_t N >
std::ostream& operator<<( std::ostream& out , const std::tr1::array< T , N > &x )
{
	if( !x.empty() ) out << x[0];
	for( size_t i=1 ; i<x.size() ; ++i ) out << "\t" << x[i];
	return out;
}

typedef boost::numeric::odeint::taylor< 3 , 10 > taylor_type;
typedef taylor_type::state_type state_type;
typedef taylor_type::derivs_type derivs_type;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;






namespace fusion = boost::fusion;

using namespace std;
using namespace boost::numeric::odeint;

using boost::numeric::odeint::taylor_adf::arg1;
using boost::numeric::odeint::taylor_adf::arg2;
using boost::numeric::odeint::taylor_adf::arg3;

int main( int argc , char **argv )
{
	cout.precision( 14 );

	taylor_type stepper;

	state_type x = {{ 10.0 , 10.0 , 10.0 }} ;

	double t = 0.0;
	double dt = 0.01;
	while( t < 10.0 )
	{
		stepper.try_step(
				fusion::make_vector
				(
						sigma * ( arg2 - arg1 ) ,
						R * arg1 - arg2 - arg1 * arg3 ,
						arg1 * arg2 - b * arg3
				) ,
				x , t , dt );

		cout << t << "\t" << x << endl;
	}

	return 0;
}
