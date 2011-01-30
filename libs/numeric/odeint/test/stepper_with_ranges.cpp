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

#include <boost/range.hpp>

#include <boost/numeric/odeint/stepper/explicit_euler.hpp>

typedef std::vector< double > state_type;
typedef std::tr1::array< double , 3 > state_type2;

struct lorenz
{
	template< class State , class Deriv >
	void operator()( const State &x_ , Deriv &dxdt_ , double t )
	{
		typename boost::range_iterator< const State >::type x = boost::begin( x_ );
		typename boost::range_iterator< Deriv >::type dxdt = boost::begin( dxdt_ );

		dxdt[0] = x[0];
		dxdt[1] = 2.0;
		dxdt[2] = 3.0;
	}
};

BOOST_AUTO_TEST_SUITE( stepper_with_ranges )

BOOST_AUTO_TEST_CASE( explicit_euler_with_range )
{
	std::vector< double > x( 3 * 2 );
	x[0] = 1.0;
	x[1] = 1.0;
	x[2] = 1.0;
	boost::numeric::odeint::explicit_euler< state_type > euler;
	euler.do_step( lorenz() , std::make_pair( x.begin() , x.begin() + 3 ) , 0.1 , 0.1 );
	BOOST_CHECK_CLOSE( x[0] , 1.1 , 1.0e-10 );
	BOOST_CHECK_CLOSE( x[1] , 1.2 , 1.0e-10 );
	BOOST_CHECK_CLOSE( x[2] , 1.3 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( explicit_euler_with_array )
{
	state_type2 x;
	boost::numeric::odeint::explicit_euler< state_type > euler;
	euler.do_step( lorenz() , x , 0.0 , 0.1 );
	euler.do_step( lorenz() , x , 0.1 , 0.1 );
}


BOOST_AUTO_TEST_SUITE_END()
