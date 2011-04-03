/*
 * main.cpp
 *
 *  Created on: Apr 2, 2011
 *      Author: karsten
 */

#include <iostream>

#include "taylor.hpp"

#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>
#include <boost/fusion/include/make_vector.hpp>

template< typename T , size_t N >
std::ostream& operator<<( std::ostream& out , const std::tr1::array< T , N > &x )
{
	if( !x.empty() ) out << x[0];
	for( size_t i=1 ; i<x.size() ; ++i ) out << "\t" << x[i];
	return out;
}

typedef boost::numeric::odeint::taylor< 3 , 11 > taylor_type;
typedef taylor_type::state_type state_type;
typedef boost::numeric::odeint::explicit_rk4< state_type > rk4_type;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}





namespace fusion = boost::fusion;

using namespace std;
using namespace boost::numeric::odeint;

using boost::numeric::odeint::taylor_adf::arg1;
using boost::numeric::odeint::taylor_adf::arg2;
using boost::numeric::odeint::taylor_adf::arg3;

int main( int argc , char **argv )
{
	taylor_type tay;
	rk4_type rk4;


	state_type x1 = {{ 10.0 , 10.0 , 10.0 }} , x2 = x1;

	const double dt = 0.1;
	double t = 0.0;
	for( size_t i=0 ; i<10000 ; ++i , t += dt )
	{
		tay.do_step(
				fusion::make_vector
				(
						sigma * ( arg2 - arg1 ) ,
						R * arg1 - arg2 - arg1 * arg3 ,
						arg1 * arg2 - b * arg3
				) ,
				x1 , t , dt );

		rk4.do_step( lorenz , x2 , t , dt );

		cout << t << "\t" << x1 << "\t" << x2 << "\n";
	}


	return 0;
}
