/*
 * main.cpp
 *
 *  Created on: Apr 2, 2011
 *      Author: karsten
 */

#include <iostream>
#include <fstream>

#include "taylor.hpp"
#include "taylor_controller.hpp"

#include <boost/numeric/odeint/stepper/explicit_error_rk54_ck.hpp>
#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>

#include <boost/fusion/include/make_vector.hpp>
#include <boost/spirit/include/phoenix.hpp>

template< typename T , size_t N >
std::ostream& operator<<( std::ostream& out , const std::tr1::array< T , N > &x )
{
	if( !x.empty() ) out << x[0];
	for( size_t i=1 ; i<x.size() ; ++i ) out << "\t" << x[i];
	return out;
}

typedef boost::numeric::odeint::taylor< 3 , 20 > taylor_type;
typedef taylor_type::state_type state_type;
typedef boost::numeric::odeint::explicit_error_rk54_ck< state_type > rk54_type;

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
namespace phoenix = boost::phoenix;

using namespace std;
using namespace boost::numeric::odeint;

using boost::numeric::odeint::taylor_adf::arg1;
using boost::numeric::odeint::taylor_adf::arg2;
using boost::numeric::odeint::taylor_adf::arg3;

struct streaming_observer
{
	std::ostream& m_out;

	streaming_observer( std::ostream &out ) : m_out( out ) { }

	template< class State >
	void operator()( const State &x , double t ) const
	{
		m_out << t;
		for( size_t i=0 ; i<x.size() ; ++i ) m_out << "\t" << x[i];
		m_out << "\n";
	}
};


int main( int argc , char **argv )
{
	double eps_abs = 1.0e-6;
	double eps_rel = 1.0e-6;
	if( argc > 2 )
	{
		eps_abs = atof( argv[1] );
		eps_rel = atof( argv[2] );
	}


	rk54_type rk54_plain;
	controlled_error_stepper< rk54_type > rk54( rk54_plain , default_error_checker< double >( eps_abs , eps_rel ) );
	taylor_controller< taylor_type > taylor_controller( eps_abs , eps_rel );

	state_type x1 = {{ 10.0 , 10.0 , 10.0 }} , x2 = x1;

	ofstream fout( "lorenz_rk54.dat" );
	size_t steps_rk54 = integrate_adaptive( rk54 , lorenz , x1 , 0.0 , 50.0 , 0.1 , streaming_observer( fout ) );
	clog << "Steps RK 54 : " << steps_rk54  << endl;;

	ofstream fout2( "lorenz_taylor.dat" );
	size_t steps_taylor = integrate_adaptive( taylor_controller ,
			fusion::make_vector
			(
					sigma * ( arg2 - arg1 ) ,
					R * arg1 - arg2 - arg1 * arg3 ,
					arg1 * arg2 - b * arg3
			) , x2 , 0.0 , 50.0 , 0.1 , streaming_observer( fout2 ) );
	clog << "Steps Taylor : " << steps_taylor << endl;

	return 0;
}
