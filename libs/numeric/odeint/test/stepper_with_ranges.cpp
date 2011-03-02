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
#include <boost/ref.hpp>

#include <boost/numeric/odeint/stepper/explicit_euler.hpp>
#include <boost/numeric/odeint/stepper/explicit_error_rk54_ck.hpp>
#include <boost/numeric/odeint/stepper/explicit_error_dopri5.hpp>
#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>
#include <boost/numeric/odeint/stepper/symplectic_euler.hpp>
#include <boost/numeric/odeint/stepper/dense_output_explicit.hpp>
#include <boost/numeric/odeint/stepper/dense_output_controlled_explicit_fsal.hpp>

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


/*
 * Useful for Hamiltonian systems
 */
struct ham_sys
{
	template< class State , class Deriv >
	void operator()( const State &x_ , Deriv &dxdt_ )
	{
		typename boost::range_iterator< Deriv >::type dxdt = boost::begin( dxdt_ );
		dxdt[0] = 1.0;
		dxdt[1] = 2.0;
		dxdt[2] = 3.0;
	}

	template< class State , class Deriv >
	void operator()( const State &x_ , const Deriv &dxdt_ )
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
	std::tr1::array< double , dim > q;
	std::tr1::array< double , dim > p;
	state_type err;

	vector_fixture( void )
//	: in( dim )
	: in() , err( 3 )
	{
		for( size_t i=0 ; i<dim ; ++i )
		{
			in[ i ] = q[i] = p[i] = double( i );
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

BOOST_AUTO_TEST_CASE( controlled_error_stepper_rk54 )
{
	double t = 0.0 , dt = 0.1;
	vector_fixture f;
	boost::numeric::odeint::controlled_error_stepper< boost::numeric::odeint::explicit_error_rk54_ck< state_type > > stepper;
	stepper.try_step( system2() , std::make_pair( f.in.begin() + 1 , f.in.begin() + 4 ) , t , dt );
	CHECK_VALUES( f.in , 0.0 , 1.1 , 2.2 , 3.3 , 4.0 , 5.0 );
}

BOOST_AUTO_TEST_CASE( controlled_error_stepper_dopri5 )
{
	double t = 0.0 , dt = 0.1;
	vector_fixture f;
	boost::numeric::odeint::controlled_error_stepper< boost::numeric::odeint::explicit_error_dopri5< state_type > > stepper;
	stepper.try_step( system2() , std::make_pair( f.in.begin() + 1 , f.in.begin() + 4 ) , t , dt );
	CHECK_VALUES( f.in , 0.0 , 1.1 , 2.2 , 3.3 , 4.0 , 5.0 );
}

BOOST_AUTO_TEST_CASE( symplectic_euler_coor_func )
{
	vector_fixture f;
	boost::numeric::odeint::symplectic_euler< state_type > euler;
	euler.do_step( ham_sys() ,
		std::make_pair( f.q.begin() + 1 , f.q.begin() + 4 ) ,
		std::make_pair( f.p.begin() + 3 , f.p.begin() + 6 ) ,
		0.0 , 0.1 );
	CHECK_VALUES( f.q , 0.0 , 1.3 , 2.4 , 3.5 , 4.0 , 5.0 );
	CHECK_VALUES( f.p , 0.0 , 1.0 , 2.0 , 3.1 , 4.2 , 5.3 );
}

BOOST_AUTO_TEST_CASE( symplectic_euler_coor_and_mom_func )
{
	vector_fixture f;
	boost::numeric::odeint::symplectic_euler< state_type > euler;
	euler.do_step( std::make_pair( ham_sys() , ham_sys() ) ,
		std::make_pair( f.q.begin() + 1 , f.q.begin() + 4 ) ,
		std::make_pair( f.p.begin() + 3 , f.p.begin() + 6 ) ,
		0.0 , 0.1 );
	CHECK_VALUES( f.q , 0.0 , 1.1 , 2.2 , 3.3 , 4.0 , 5.0 );
	CHECK_VALUES( f.p , 0.0 , 1.0 , 2.0 , 3.1 , 4.2 , 5.3 );
}

BOOST_AUTO_TEST_CASE( dense_output_euler_with_ranges )
{
	using namespace boost::numeric::odeint;
	vector_fixture f;
	dense_output_explicit< explicit_euler< state_type > > stepper;
	stepper.initialize( std::make_pair( f.in.begin() + 1, f.in.begin() + 4 ) , 0.0 , 0.1 );
	stepper.do_step( system1() );
	stepper.calc_state( 0.05 , std::make_pair( f.in.begin() + 1 ,f.in.begin() +4 ) );
	CHECK_VALUES( f.in , 0.0 , 1.05 , 2.1 , 3.15 , 4.0 , 5.0 );
}

BOOST_AUTO_TEST_CASE( dense_output_dopri5_with_ranges )
{
	using namespace boost::numeric::odeint;
	vector_fixture f;
	dense_output_controlled_explicit_fsal<
		controlled_error_stepper<
			explicit_error_dopri5< state_type >
		> > stepper;
	stepper.initialize( std::make_pair( f.in.begin() + 1, f.in.begin() + 4 ) , 0.0 , 0.1 );
	stepper.do_step( system2() );
	stepper.calc_state( 0.05 , std::make_pair( f.in.begin() + 1 ,f.in.begin() +4 ) );
	CHECK_VALUES( f.in , 0.0 , 1.05 , 2.1 , 3.15 , 4.0 , 5.0 );
}




BOOST_AUTO_TEST_SUITE_END()
