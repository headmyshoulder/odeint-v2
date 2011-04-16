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
	cout.precision( 14 );

	taylor_type t;

	state_type in = {{ 10.0 , 10.0 , 10.0 }} , dxdt = {{ 0.0 , 0.0 , 0.0 }} , xerr = {{ 0.0 , 0.0 , 0.0 }} , out = {{ 0.0 ,0.0 , 0.0 }};

	lorenz( in , dxdt , 0.0 );

	cout << in << endl;
	cout << dxdt << endl << endl;

	t.do_step(
			fusion::make_vector
			(
					sigma * ( arg2 - arg1 ) ,
					R * arg1 - arg2 - arg1 * arg3 ,
					arg1 * arg2 - b * arg3
			) ,
			in , 0.0 , out , 0.1 , xerr );



	cout << endl;
	cout << in << endl;
	cout << xerr << endl;
	cout << out << endl << endl;
	const derivs_type &derivs = t.get_last_derivs();
	for( size_t i=0 ; i<derivs.size() ; ++i )
		cout << derivs[i] << endl;

	return 0;
}
