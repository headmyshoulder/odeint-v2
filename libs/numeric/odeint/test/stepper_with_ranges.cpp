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
#include <boost/numeric/odeint/stepper/explicit_error_rk54_ck.hpp>
#include <boost/numeric/odeint/stepper/explicit_error_dopri5.hpp>

typedef std::vector< double > state_type;
typedef std::tr1::array< double , 3 > state_type2;


/*
 * The two systems are needed, since for steppers with more than
 * one internal step it is difficult to calculate the exact result
 *
 * system1 is suited for euler
 */
struct system1
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

/*
 * system2 is suited for all steppers, it allows you to calculate the result analytically.
 */
struct system2
{
	template< class State , class Deriv >
	void operator()( const State &x_ , Deriv &dxdt_ , double t )
	{
		typename boost::range_iterator< Deriv >::type dxdt = boost::begin( dxdt_ );

		dxdt[0] = 1.0;
		dxdt[1] = 2.0;
		dxdt[2] = 3.0;
	}

	template< class State , class Deriv >
	void operator()( const State &x_ , const Deriv &dxdt_ , double t )
	{
		typename boost::range_iterator< Deriv >::type dxdt = boost::begin( dxdt_ );

		dxdt[0] = 1.0;
		dxdt[1] = 2.0;
		dxdt[2] = 3.0;
	}
};


struct vector_fixture
{
	const static size_t dim = 6;
	std::tr1::array< double , dim > in;
	state_type err;

	vector_fixture( void )
//	: in( dim )
	: in() , err( 3 )
	{
		for( size_t i=0 ; i<dim ; ++i )
		{
			in[ i ] = double( i );
		}
		for( size_t i=0 ; i<3 ; ++i )
		{
			err[i] = double( i ) * 10.0;
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



BOOST_AUTO_TEST_SUITE( stepper_with_ranges )

BOOST_AUTO_TEST_CASE( explicit_euler_with_range_v1 )
{
	vector_fixture f;
	boost::numeric::odeint::explicit_euler< state_type > euler;
	euler.do_step( system1() , std::make_pair( f.in.begin() + 1 , f.in.begin() + 4 ) , 0.1 , 0.1 );
	CHECK_VALUES( f.in , 0.0 , 1.1 , 2.2 , 3.3 , 4.0 , 5.0 );
}

BOOST_AUTO_TEST_CASE( explicit_error_k54_with_range_v1 )
{
	vector_fixture f;
	boost::numeric::odeint::explicit_error_rk54_ck< state_type > rk54;
	rk54.do_step( system2() , std::make_pair( f.in.begin() + 1 , f.in.begin() + 4 ) , 0.1 , 0.1 );
	CHECK_VALUES( f.in , 0.0 , 1.1 , 2.2 , 3.3 , 4.0 , 5.0 );
}

BOOST_AUTO_TEST_CASE( explicit_error_k54_with_range_v5 )
{
	vector_fixture f;
	boost::numeric::odeint::explicit_error_rk54_ck< state_type > rk54;
	rk54.do_step( system2() , std::make_pair( f.in.begin() + 1 , f.in.begin() + 4 ) , 0.1 , 0.1 , f.err );
	CHECK_VALUES( f.in , 0.0 , 1.1 , 2.2 , 3.3 , 4.0 , 5.0 );
}


BOOST_AUTO_TEST_CASE( explicit_error_dopri5_with_range_v1 )
{
	vector_fixture f;
	boost::numeric::odeint::explicit_error_dopri5< state_type > dopri5;
	dopri5.do_step( system2() , std::make_pair( f.in.begin() + 1 , f.in.begin() + 4 ) , 0.1 , 0.1 );
	CHECK_VALUES( f.in , 0.0 , 1.1 , 2.2 , 3.3 , 4.0 , 5.0 );
}

BOOST_AUTO_TEST_CASE( explicit_error_dopri5_with_range_v5 )
{
	vector_fixture f;
	boost::numeric::odeint::explicit_error_dopri5< state_type > dopri5;
	dopri5.do_step( system2() , std::make_pair( f.in.begin() + 1 , f.in.begin() + 4 ) , 0.1 , 0.1 , f.err );
	CHECK_VALUES( f.in , 0.0 , 1.1 , 2.2 , 3.3 , 4.0 , 5.0 );
}



BOOST_AUTO_TEST_SUITE_END()
