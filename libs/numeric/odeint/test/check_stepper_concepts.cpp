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
#include <boost/numeric/odeint/algebra/tr1_array_resize.hpp>

#include "vector_space_1d.hpp"
#include "gsl_vector_adaptor.hpp"

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


//template< class ErrorStepper >
//void check_error_stepper_concept(
//    ErrorStepper &stepper ,
//    typename ErrorStepper::order_type order_error_step ,
//    typename ErrorStepper::order_type order_error )
//{
//    typedef ErrorStepper stepper_type;
//    typedef typename stepper_type::container_type container_type;
//    typedef typename stepper_type::traits_type traits_type;
//    typedef typename stepper_type::value_type value_type;
//    typedef typename stepper_type::order_type order_type;
//    typedef typename stepper_type::time_type time_type;
//
//    constant_system< container_type > con;
//
//    BOOST_CHECK_EQUAL( order_error_step , stepper.order_error_step() );
//    BOOST_CHECK_EQUAL( order_error , stepper.order_error() );
//
//    container_type x( 1 , 0.0 ) , xerr( 1 , 0.0 );
//    stepper.adjust_size( x );
//
//    stepper.do_step( con , x , 0.0 , 0.1 , xerr );
//    BOOST_CHECK_SMALL( fabs( x[0] - 0.1 ) , eps );
//    BOOST_CHECK_SMALL( fabs( xerr[0] ) , eps );
//
//    container_type dxdt( 1 , 1.0 );
//    stepper.do_step( con , x , dxdt , 0.0 , 0.1 , xerr );
//    BOOST_CHECK_SMALL( fabs( x[0] - 0.2 ) , eps );
//    BOOST_CHECK_SMALL( fabs( xerr[0] ) , eps );
//
//    stepper_type stepper2( x );
//    stepper_type stepper3;
//}
//
//
//
//
//
//template< class ControlledErrorStepper >
//void check_controlled_stepper_concept(
//    ControlledErrorStepper &stepper
//    )
//{
//    typedef ControlledErrorStepper stepper_type;
//    typedef typename stepper_type::container_type container_type;
//    typedef typename stepper_type::traits_type traits_type;
//    typedef typename stepper_type::value_type value_type;
//    typedef typename stepper_type::order_type order_type;
//    typedef typename stepper_type::time_type time_type;
//
////    constant_system< container_type > con;
//
//    container_type x( 1 , 0.0 );
//    stepper.adjust_size( x );
//}
//
//
//
//
//
//void test_euler_concept()
//{
//    stepper_euler< std::vector<double> > stepper;
//    check_stepper_concept( stepper , 1 );
//}
//
//
//
//void test_half_step_euler_concept()
//{
//    stepper_half_step< stepper_euler< std::vector< double > > > stepper;
//    check_stepper_concept( stepper , 1 );
//    check_error_stepper_concept( stepper , 1 , 2 );
//}
//
///*
//void test_midpoint_concept()
//{
//    stepper_midpoint< std::vector< double > > stepper;
//    stepper.set_step_number( 4 );
//    unsigned short step_number = stepper.get_step_number();
//    step_number = 5; // no warnings
//    check_stepper_concept( stepper , 2 );
//}
//
//void test_rk4_classical_concept()
//{
//    stepper_rk4_classical< std::vector<double> > stepper;
//    check_stepper_concept( stepper , 4 );
//}
//
//void test_rk4_concept()
//{
//    stepper_rk4< std::vector<double> > stepper;
//    check_stepper_concept( stepper , 4 );
//}
//
//void test_rk5_ck_concept()
//{
//    stepper_rk5_ck< std::vector<double> > stepper;
//    check_error_stepper_concept( stepper , 5 , 5 );
//}
//
//void test_rk78_fehlberg_concept()
//{
//    stepper_rk78_fehlberg< std::vector<double> > stepper;
//    check_stepper_concept( stepper , 8 );
//    check_error_stepper_concept( stepper , 7 , 8 );
//}
//*/
//
//void test_controlled_stepper_standard_concept()
//{
//    typedef stepper_euler< std::vector< double > > stepper_type;
//    typedef controlled_stepper_standard< stepper_type > controlled_stepper_type;
//
//    controlled_stepper_type stepper( 1.0 , 1.0 , 1.0 , 1.0 );
//    check_controlled_stepper_concept( stepper );
//    }
//


void test_euler_with_vector( void )
{
	state_type1 x( 1 , 0.0 );
	explicit_euler< state_type1 > euler;
	check_stepper_concept( euler , constant_system1 , x );
}

void test_euler_with_array( void )
{
	state_type4 x;
	x[0] = 0.0;
	explicit_euler< state_type4 > euler;
	check_stepper_concept( euler , constant_system4 , x );
}

void test_runge_kutta_error_ck_with_vector( void )
{
	state_type1 x( 1 , 0.0 );
	state_type1 xerr( 1 , 0.0 );
	runge_kutta_error_ck< state_type1 > rk_ck;
	check_error_stepper_concept( rk_ck , constant_system1 , x , xerr );
}

test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE("check stepper concepts");



    test->add( BOOST_TEST_CASE( &test_euler_with_vector ) );

//    test->add( BOOST_TEST_CASE( &test_euler_concept ) );

    return test;
}
