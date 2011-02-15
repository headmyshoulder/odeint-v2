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

	template< class State , class Deriv >
	void operator()( const State &x_ , const Deriv &dxdt_ , double t )
	{
		typename boost::range_iterator< const State >::type x = boost::begin( x_ );
		typename boost::range_iterator< Deriv >::type dxdt = boost::begin( dxdt_ );

		dxdt[0] = x[0];
		dxdt[1] = 2.0;
		dxdt[2] = 3.0;
	}

};

struct vector_fixture
{
	const static size_t dim = 6;
	std::vector< double > in , out , dxdt ;
	boost::numeric::odeint::explicit_euler< state_type > euler;

	vector_fixture( void )
	: in( dim ) , out( dim ) , dxdt( dim )
	{
		for( size_t i=0 ; i<dim ; ++i )
		{
			in[ i ] = double( i );
			out[ i ] = double( i ) + 10.0 ;
			dxdt[ i ] = double( i ) + 100.0 ;
		}
	}

	~vector_fixture( void )
	{
	}
};

#define CHECK_VALUES( x , x0 , x1 , x2 , x3 , x4 , x5 ) \
	BOOST_CHECK_CLOSE( x[0] , x0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[1] , x1 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[2] , x2 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[3] , x3 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[4] , x4 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[5] , x5 , 1.0e-8 )

#define CHECK_IN_DEFAULT( x ) \
	BOOST_CHECK_CLOSE( x[0] , 0.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[1] , 1.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[2] , 2.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[3] , 3.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[4] , 4.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[5] , 5.0 , 1.0e-8 )

#define CHECK_OUT_DEFAULT( x ) \
	BOOST_CHECK_CLOSE( x[0] , 10.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[1] , 11.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[2] , 12.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[3] , 13.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[4] , 14.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[5] , 15.0 , 1.0e-8 )

#define CHECK_DXDT_DEFAULT( x ) \
	BOOST_CHECK_CLOSE( x[0] , 100.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[1] , 101.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[2] , 102.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[3] , 103.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[4] , 104.0 , 1.0e-8 ); \
	BOOST_CHECK_CLOSE( x[5] , 105.0 , 1.0e-8 )



BOOST_AUTO_TEST_SUITE( stepper_with_ranges )

BOOST_AUTO_TEST_CASE( explicit_euler_with_range_v1 )
{
	vector_fixture f;
	f.euler.do_step( lorenz() , std::make_pair( f.in.begin() + 1 , f.in.begin() + 4 ) , 0.1 , 0.1 );
	CHECK_VALUES( f.in , 0.0 , 1.1 , 2.2 , 3.3 , 4.0 , 5.0 );
	CHECK_OUT_DEFAULT( f.out );
	CHECK_DXDT_DEFAULT( f.dxdt );
}

BOOST_AUTO_TEST_CASE( explicit_euler_with_range_v2 )
{
	vector_fixture f;
	lorenz()( std::make_pair( f.in.begin() , f.in.begin() + 3 ) ,
			  std::make_pair( f.dxdt.begin() + 3 , f.dxdt.begin() + 6 ) , 0.0 );
	f.euler.do_step( lorenz() , std::make_pair( f.in.begin() , f.in.begin() + 3 ) ,
								std::make_pair( f.dxdt.begin() + 3 , f.dxdt.begin() + 6 ) , 0.1 , 0.1 );
	CHECK_VALUES( f.in , 0.0 , 1.2 , 2.3 , 3.0 , 4.0 , 5.0 );
	CHECK_OUT_DEFAULT( f.out );
	CHECK_VALUES( f.dxdt , 100.0 , 101.0 , 102.0 , 0.0 , 2.0 , 3.0 );
}

BOOST_AUTO_TEST_CASE( explicit_euler_with_range_v3 )
{
	vector_fixture f;
	f.euler.do_step( lorenz() , std::make_pair( f.in.begin() + 2 , f.in.begin() + 5 ) , 0.1 ,
			std::make_pair( f.out.begin() , f.out.begin() +3 ) , 0.1 );
	CHECK_IN_DEFAULT( f.in );
	CHECK_VALUES( f.out , 2.2 , 3.2 , 4.3 , 13.0 , 14.0 , 15.0 );
	CHECK_DXDT_DEFAULT( f.dxdt );
}

BOOST_AUTO_TEST_CASE( explicit_euler_with_range_v4 )
{
	vector_fixture f;
	f.euler.do_step( lorenz() , std::make_pair( f.in.begin() + 2 , f.in.begin() + 5 ) , 0.1 ,
			std::make_pair( f.out.begin() , f.out.begin() +3 ) , 0.1 );
	CHECK_IN_DEFAULT( f.in );
	CHECK_VALUES( f.out , 2.2 , 3.2 , 4.3 , 13.0 , 14.0 , 15.0 );
	CHECK_DXDT_DEFAULT( f.dxdt );
}



BOOST_AUTO_TEST_CASE( explicit_euler_with_array )
{
	state_type2 x;
	boost::numeric::odeint::explicit_euler< state_type > euler;
	euler.do_step( lorenz() , x , 0.0 , 0.1 );
	euler.do_step( lorenz() , x , 0.1 , 0.1 );
}


BOOST_AUTO_TEST_SUITE_END()
