/* Boost stepper_euler.cpp test file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file tests the use of the euler stepper
  
 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <vector>
#include <cmath>
#include <tr1/array>

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/algebra/gsl_vector_adaptor.hpp>

#include "vector_space_1d.hpp"




using std::vector;

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

typedef std::vector< double > state_type1;
typedef gsl_vector state_type2;
typedef vector_space_1d< double > state_type3;
typedef std::tr1::array< double , 1 > state_type4;

void constant_system1( const state_type1 &x , state_type1 &dxdt , double t )
{
	dxdt[0] = 1.0;
}

void constant_system2( const state_type2 &x , state_type2 &dxdt , double t )
{
	gsl_vector_set( &dxdt , 0 , 1.0 );
}

void constant_system3( const state_type3 &x , state_type3 &dxdt , double t )
{
	dxdt.m_x = 1.0;
}

void constant_system4( const state_type4 &x , state_type4 &dxdt , double t )
{
	dxdt[0] = 1.0;
}



const double eps = 1.0e-14;


template< class Stepper , class System >
void check_stepper_concept( Stepper &stepper , System system , typename Stepper::state_type &x )
{
    typedef Stepper stepper_type;
    typedef typename stepper_type::state_type container_type;
    typedef typename stepper_type::order_type order_type;
    typedef typename stepper_type::time_type time_type;

    stepper.do_step( system , x , 0.0 , 0.1 );
    double xval = * boost::begin( x );
    BOOST_CHECK_SMALL( fabs( xval - 0.1 ) , eps );
}

template< class Stepper , class System >
void check_error_stepper_concept( Stepper &stepper , System system ,
									  typename Stepper::state_type &x , typename Stepper::state_type &xerr )
{
    typedef Stepper stepper_type;
    typedef typename stepper_type::state_type container_type;
    typedef typename stepper_type::order_type order_type;
    typedef typename stepper_type::time_type time_type;

    stepper.do_step( system , x , 0.0 , 0.1 , xerr);
    double xval = * boost::begin( x );
    BOOST_CHECK_SMALL( fabs( xval - 0.1 ) , eps );
}

void test_euler_with_vector( void )
{
	state_type1 x( 1 , 0.0 );
	explicit_euler< state_type1 > euler;
	check_stepper_concept( euler , constant_system1 , x );
}

void test_euler_with_gsl_vector( void )
{
	state_type2 *x = gsl_vector_alloc( 1 );
	explicit_euler< state_type2 > euler;
//	check_stepper_concept( euler , constant_system2 , *x );
	gsl_vector_free( x );
}

void test_euler_with_array( void )
{
	state_type4 x;
	x[0] = 0.0;
	explicit_euler< state_type4 > euler;
	check_stepper_concept( euler , constant_system4 , x );
}

//void test_runge_kutta_error_ck_with_vector( void )
//{
//	state_type1 x( 1 , 0.0 );
//	state_type1 xerr( 1 , 0.0 );
//	runge_kutta_error_ck< state_type1 > rk_ck;
//	check_error_stepper_concept( rk_ck , constant_system1 , x , xerr );
//}

test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE("check stepper concepts");



    test->add( BOOST_TEST_CASE( &test_euler_with_vector ) );
    test->add( BOOST_TEST_CASE( &test_euler_with_array ) );

//    test->add( BOOST_TEST_CASE( &test_euler_with_gsl_vector ) );



    return test;
}
