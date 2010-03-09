/* Boost stepper_euler.cpp test file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file tests the use of the euler stepper
  
 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <vector>
#include <tr1/array>

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/stepper_euler.hpp>
#include <boost/numeric/odeint/stepper_half_step.hpp>
#include <boost/numeric/odeint/stepper_midpoint.hpp>
#include <boost/numeric/odeint/stepper_rk4_classical.hpp>
#include <boost/numeric/odeint/stepper_rk4.hpp>
#include <boost/numeric/odeint/stepper_rk5_ck.hpp>
#include <boost/numeric/odeint/stepper_rk78_fehlberg.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;




template< class Container >
struct constant_system
{
    void operator()( const Container &x , Container &dxdt , double t )
    {
        dxdt[0] = 1.0;
    }
};


const double eps = 1.0e-14;

template< class Stepper >
void check_stepper_concept( Stepper &stepper ,
                            typename Stepper::order_type order_step )
{
    typedef Stepper stepper_type;
    typedef typename stepper_type::container_type container_type;
    typedef typename stepper_type::traits_type traits_type;
    typedef typename stepper_type::value_type value_type;
    typedef typename stepper_type::order_type order_type;
    typedef typename stepper_type::time_type time_type;

    constant_system< container_type > con;

    BOOST_CHECK_EQUAL( order_step , stepper.order_step() );

    container_type x( 1 , 0.0 ) ;
    stepper.adjust_size( x );
    stepper.do_step( con , x , 0.0 , 0.1 );
    BOOST_CHECK_SMALL( fabs( x[0] - 0.1 ) , eps );

    container_type dxdt( 1 , 1.0 );
    stepper.do_step( con , x , dxdt , 0.0 , 0.1 );
    BOOST_CHECK_SMALL( fabs( x[0] - 0.2 ) , eps );

    stepper_type stepper2( x );
    stepper_type stepper3;
}




template< class ErrorStepper >
void check_error_stepper_concept(
    ErrorStepper &stepper ,
    typename ErrorStepper::order_type order_error_step ,
    typename ErrorStepper::order_type order_error )
{
    typedef ErrorStepper stepper_type;
    typedef typename stepper_type::container_type container_type;
    typedef typename stepper_type::traits_type traits_type;
    typedef typename stepper_type::value_type value_type;
    typedef typename stepper_type::order_type order_type;
    typedef typename stepper_type::time_type time_type;

    constant_system< container_type > con;

    BOOST_CHECK_EQUAL( order_error_step , stepper.order_error_step() );
    BOOST_CHECK_EQUAL( order_error , stepper.order_error() );

    container_type x( 1 , 0.0 ) , xerr( 1 , 0.0 );
    stepper.adjust_size( x );

    stepper.do_step( con , x , 0.0 , 0.1 , xerr );
    BOOST_CHECK_SMALL( fabs( x[0] - 0.1 ) , eps );
    BOOST_CHECK_SMALL( fabs( xerr[0] ) , eps );

    container_type dxdt( 1 , 1.0 );
    stepper.do_step( con , x , dxdt , 0.0 , 0.1 , xerr );
    BOOST_CHECK_SMALL( fabs( x[0] - 0.2 ) , eps );
    BOOST_CHECK_SMALL( fabs( xerr[0] ) , eps );

    stepper_type stepper2( x );
    stepper_type stepper3;
}





void test_euler_concept()
{
    stepper_euler< std::vector<double> > stepper;
    check_stepper_concept( stepper , 1 );
}



void test_half_step_euler_concept()
{
    stepper_half_step< stepper_euler< std::vector< double > > > stepper;
    check_stepper_concept( stepper , 1 );
    check_error_stepper_concept( stepper , 1 , 2 );
}

void test_midpoint_concept()
{
    stepper_midpoint< std::vector< double > > stepper;
    stepper.set_step_number( 4 );
    unsigned short step_number = stepper.get_step_number();
    step_number = 5; // no warnings
    check_stepper_concept( stepper , 2 );
}

void test_rk4_classical_concept()
{
    stepper_rk4_classical< std::vector<double> > stepper;
    check_stepper_concept( stepper , 4 );
}

void test_rk4_concept()
{
    stepper_rk4< std::vector<double> > stepper;
    check_stepper_concept( stepper , 4 );
}

void test_rk5_ck_concept()
{
    stepper_rk5_ck< std::vector<double> > stepper;
    check_error_stepper_concept( stepper , 5 , 5 );
}

void test_rk78_fehlberg_concept()
{
    stepper_rk78_fehlberg< std::vector<double> > stepper;
    check_stepper_concept( stepper , 8 );
//    check_error_stepper_concept( stepper , 7 , 8 );
}




test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE("check stepper concepts");

    test->add( BOOST_TEST_CASE( &test_euler_concept ) );
    test->add( BOOST_TEST_CASE( &test_half_step_euler_concept ) );
    test->add( BOOST_TEST_CASE( &test_midpoint_concept ) );
    test->add( BOOST_TEST_CASE( &test_rk4_classical_concept ) );
    test->add( BOOST_TEST_CASE( &test_rk4_concept ) );
    test->add( BOOST_TEST_CASE( &test_rk5_ck_concept ) );
    test->add( BOOST_TEST_CASE( &test_rk78_fehlberg_concept ) );

    return test;
}
