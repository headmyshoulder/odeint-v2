/*
 * adams_bashforth.cpp
 *
 *  Created on: May 18, 2011
 *      Author: karsten
 */

#include <ostream>
#include <iostream>

#include <boost/array.hpp>
#include <boost/numeric/odeint/stepper/adams_bashforth_moulton.hpp>
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef double value_type;
typedef boost::array< value_type , 3 > state_type;

struct lorenz
{
	template< class State , class Deriv , class Value >
	void operator()( const State &_x , Deriv &_dxdt , const Value &dt ) const
	{
        const value_type sigma = 10.0;
        const value_type R = 28.0;
        const value_type b = 8.0 / 3.0;

		typename boost::range_iterator< const State >::type x = boost::begin( _x );
		typename boost::range_iterator< Deriv >::type dxdt = boost::begin( _dxdt );

        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = x[0]*x[1] - b * x[2];
	}
};

struct writing_observer
{
	std::ostream &m_out;
	writing_observer( std::ostream &out ) : m_out( out ) { }

	template< class State , class Value >
	void operator()( const State &x , const Value &t )
	{
		m_out << t << "\t" << x[0] << "\t" << x[1] << "\t" << x[2] << "\n";
	}
};


int main( int argc , char **argv )
{
	state_type x = {{ 10.0 , 10.0 , 10.0 }};

	const double dt = 0.01;
	double t = 0.0;

	adams_bashforth_moulton< 5 , state_type > stepper;
	stepper.initialize( lorenz() , x , t , dt );

	integrate_n_steps( stepper , lorenz() , x , t , dt , 10000 , writing_observer( cout ) );
}
