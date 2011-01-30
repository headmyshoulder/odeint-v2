/*
 * stepper_with_ranges.cpp
 *
 *  Created on: Jan 30, 2011
 *      Author: karsten
 */

#define BOOST_TEST_MODULE odeint_stepper_with_ranges

#include <boost/test/unit_test.hpp>

#include <vector>
#include <tr1/array>
#include <utility>

#include <boost/numeric/odeint/stepper/explicit_euler.hpp>

typedef std::vector< double > state_type;
typedef std::tr1::array< double , 3 > state_type2;

struct lorenz
{
	template< class State , class Deriv >
	void operator()( const State &x , Deriv &dxdt , double t )
	{
		const double sigma = 10.0;
		const double R = 28.0;
		const double b = 8.0 / 3.0;

//		dxdt[0] = sigma * ( x[1] - x[0] );
//		dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
//		dxdt[2] = x[0]*x[1] - b * x[2];
	}
};

BOOST_AUTO_TEST_SUITE( stepper_with_ranges )

BOOST_AUTO_TEST_CASE( explicit_euler_with_range )
{
	std::vector< double > x( 3 * 2 );
	boost::numeric::odeint::explicit_euler< state_type > euler;
	std::pair< std::vector< double >::iterator , std::vector< double >::iterator > r( x.begin() , x.begin() + 3 );
	euler.do_step( lorenz() , r , 0.0 , 0.1 );
	euler.do_step( lorenz() , std::make_pair( x.begin() + 3 , x.begin() + 6 ) , 0.1 , 0.1 );
}

BOOST_AUTO_TEST_CASE( explicit_euler_with_array )
{
	state_type2 x;
	boost::numeric::odeint::explicit_euler< state_type > euler;
	euler.do_step( lorenz() , x , 0.0 , 0.1 );
	euler.do_step( lorenz() , x , 0.1 , 0.1 );
}


BOOST_AUTO_TEST_SUITE_END()
