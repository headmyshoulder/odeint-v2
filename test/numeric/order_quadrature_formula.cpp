/* Boost numeric test for orders of quadrature formulas test file

 Copyright 2012 Mario Mulansky
 Copyright 2012 Karsten Ahnert

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

// disable checked iterator warning for msvc
#include <boost/config.hpp>
/*
#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif
*/
#define BOOST_TEST_MODULE order_quadrature_formula

#include <iostream>
#include <cmath>

#include <boost/test/unit_test.hpp>

#include <boost/mpl/vector.hpp>

#include <boost/numeric/odeint.hpp>

#include <boost/numeric/ublas/vector.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;
namespace mpl = boost::mpl;

typedef double value_type;
typedef value_type time_type;
typedef value_type state_type;

BOOST_AUTO_TEST_SUITE( order_of_convergence_test )

value_type tolerance = 1.0e-13;

struct monome
{
    int exponent;
    monome() : exponent( 0 ){};
    void operator()( const state_type &x , state_type &dxdt , const time_type t ){
	dxdt = ( 1.0 + exponent )*pow( t, exponent );
    }
};

monome rhs;

/* generic test for all steppers that support integrate_const */
template< class Stepper >
struct integrate_const_test
{
    double maxError;
    int estimatedOrder;
    void operator()( int nSteps = 1 )
    {
	estimateOrder();
	std::cout << boost::format( "%-20i%-20i%-20E\n" )
	    % estimatedOrder %  definedOrder() % maxError;

	BOOST_REQUIRE_EQUAL( estimatedOrder, definedOrder() );
    }

    const int definedOrder(){
	const Stepper stepper;
	return stepper.order();
    }
    /*
    the order of the stepper is estimated by trying to solve the ODE
    x'(t) = t^p
    until the errors are too big to be justified by finite precision.
    the first value p for which the problem is *not* solved with
    precision `tolerance` is the estimate for the order of the scheme.
     */
    void estimateOrder( int nSteps = 1 )
    {
	state_type x;
	double t;
	maxError = 0;
	rhs.exponent = -1;
	do{
	    // begin with x'(t) = ( t^0 )/1 = 1
	    //         => x (t) =   t
	    // then use   x'(t) = ( t^1 )/2 = t/2
	    //         => x (t) =   t^2
	    // ...
	    rhs.exponent++;
	    Stepper stepper;
	    x = 0.0;
	    t = 0;
	    for ( int i = 0; i < nSteps; i++ ){
		stepper.do_step( rhs, x, t, 1.0/nSteps );
		t += 1.0/nSteps;
		// compute the error using the exact solution x(t)=t^(1+p)
		value_type error = fabs( x - pow( t, ( 1.0 + rhs.exponent ) ) );
		maxError = fmax( error, maxError );
	    }
	}
	while ( maxError < tolerance );
	// return the first exponent for which the test failed
	estimatedOrder = rhs.exponent;
    }
};


template< class Stepper >
struct integrate_const_test_initialize
{
    int estimatedOrder;
    value_type maxError;

    void operator()( int nSteps = 16 )
    {
	estimateOrder( nSteps );
	std::cout << boost::format( "%-20i%-20i%-30E\n" )
	    % estimatedOrder %  definedOrder() %  maxError;

	BOOST_REQUIRE_EQUAL( estimatedOrder, definedOrder() );
    }

    const int definedOrder(){
	const Stepper stepper;
	return stepper.order();
    }
    /*
    just like the other version of estimateOrder(), but with
    initization performed by the fehlberg stepper ( order 8 )

    if the default initialization ( runge kutta 4 ) is used,
    the estimated order will never be greater than 4.
    */
    void estimateOrder( int nSteps = 16 )
    {
	state_type x;
	time_type t;
	const time_type dt = 1.0/nSteps;
	maxError = 0.0;
	rhs.exponent = -1;
	do{
	    rhs.exponent++;
	    // construct the stepper inside the for loop to reset the
	    // step storage
	    Stepper stepper;
	    x = 0.0;
	    t = 0.0;
	    // use a high order method to initialize.
	    stepper.initialize( runge_kutta_fehlberg78< state_type >(),
				rhs, x, t, dt );
	    while ( t < 1 ){
		stepper.do_step( rhs, x, t, 1.0/nSteps );
		t += 1.0/nSteps;
		// compute the error using the exact solution x(t)=t^(1+p)
		value_type error = fabs( x - pow( t, 1.0 + rhs.exponent) );
		maxError = fmax( error, maxError );
	    }
	} while( maxError < tolerance );
	// return the first exponent for which the test failed
	estimatedOrder = rhs.exponent;
    }
};


typedef mpl::vector<
    euler< state_type > ,
    modified_midpoint< state_type > ,
    runge_kutta4< state_type > ,
    runge_kutta4_classic< state_type > ,
    runge_kutta_cash_karp54_classic< state_type > ,
    runge_kutta_cash_karp54< state_type > ,
    runge_kutta_dopri5< state_type > ,
    runge_kutta_fehlberg78< state_type >
    > runge_kutta_steppers;

typedef mpl::vector<
    adams_bashforth< 2, state_type > ,
    adams_bashforth< 3, state_type > ,
    adams_bashforth< 4, state_type > ,
    adams_bashforth< 5, state_type > ,
    adams_bashforth< 6, state_type > ,
    adams_bashforth< 7, state_type > ,
    adams_bashforth< 8, state_type > ,
    adams_bashforth_moulton< 2, state_type > ,
    adams_bashforth_moulton< 3, state_type > ,
    adams_bashforth_moulton< 4, state_type > ,
    adams_bashforth_moulton< 5, state_type > ,
    adams_bashforth_moulton< 6, state_type > ,
    adams_bashforth_moulton< 7, state_type > ,
    adams_bashforth_moulton< 8, state_type >
    > adams_steppers;

BOOST_AUTO_TEST_CASE( print_header )
{
    std::cout << boost::format( "%-20s%-20s%-30s\n" )
	% "Estimated order" % "defined order"
	% "maxError";
}

BOOST_AUTO_TEST_CASE_TEMPLATE( runge_kutta_test , Stepper, runge_kutta_steppers )
{
    integrate_const_test< Stepper > tester;
    tester();
}

BOOST_AUTO_TEST_CASE( print_seperator )
{
    std::cout << "-------------------"
	      << " ABM STEPPERS "
	      << "-------------------"
	      << std::endl;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( adams_moultion_test , Stepper, adams_steppers )
{
    integrate_const_test_initialize< Stepper > tester;
    tester();
}

BOOST_AUTO_TEST_SUITE_END()
