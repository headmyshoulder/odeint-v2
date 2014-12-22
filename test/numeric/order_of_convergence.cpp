/* Boost numeric test for orders of convergence file

 Copyright 2012 Mario Mulansky
 Copyright 2012 Karsten Ahnert

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

// disable checked iterator warning for msvc
#include <boost/config.hpp>

#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif

#define BOOST_TEST_MODULE order_of_convergence

#include <iostream>
#include <cmath>
#include <string>

#include <boost/array.hpp>

#include <boost/test/unit_test.hpp>

#include <boost/mpl/vector.hpp>

#include <boost/numeric/odeint.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;
namespace mpl = boost::mpl;

typedef double value_type;
typedef boost::array< double , 1 > state_type;

BOOST_AUTO_TEST_SUITE( order_of_convergence_test )

int p;
double tolerance = 1.0e-13;

/* generic test for all steppers that support integrate_const */
template< class Stepper >
struct integrate_const_test
{
    double error;
    void operator()( int nSteps = 1 )
    {
	Stepper stepper;
	std::cout << boost::format( "%-20i%-20i%-20E\n" )
	    % estimatedOrder() %  definedOrder() % error;
	BOOST_REQUIRE( estimatedOrder() == definedOrder() );
    }

    const int definedOrder(){
	const Stepper stepper;
	return stepper.order();
    }

    const int estimatedOrder( int nSteps = 1 )
    {
	state_type x;
	double t;
	for ( p = 0; p < 20; p++ )
	    {
		Stepper stepper;
		x[0] = 1.0;
		t = 0;
		for ( int i = 0; i < nSteps; i++ ){
		    stepper.do_step( rhs, x, t, 1.0/nSteps );
		    t += 1.0/nSteps;
		    error = fabs( x[0] - pow( t, (1.0+p) )/(1.0 + p) - 1.0 );
		    if ( error >tolerance )
			return p;
		}
	    }
	return p;
    }


private:
    static void rhs( const state_type &x , state_type &dxdt , const double t )
    {
	dxdt[0] = pow( t, p );
    }
};


template< class Stepper >
struct integrate_const_test_initialize
{
    double error;
    void operator()( int nSteps = 16 )
    {
	std::cout << boost::format( "%-20i%-20i%-20E" ) % estimatedOrder()
	    %  definedOrder() % error << std::endl;
	BOOST_REQUIRE( estimatedOrder() == definedOrder() );
    }

    const int definedOrder(){
	const Stepper stepper;
	return stepper.order();
    }

    const int estimatedOrder( int nSteps = 16 )
    {
	state_type x;
	double t;
	for ( p = 0; p < 10; p++ )
	    {
		Stepper stepper;
		x[0] = 1.0;
		t = 0;
		// use a high order method to initialize.
		stepper.initialize( runge_kutta_fehlberg78< state_type >(),
					 rhs, x, t, 1.0/nSteps );
		while ( t < 1 ){
		    stepper.do_step( rhs, x, t, 1.0/nSteps );
		    t += 1.0/nSteps;
		    error = fabs( x[0]-pow( t, 1.0 + p)/(1.0+p) - 1.0 );
		    if (error > tolerance )
			return p;
		}
	    }	
	return p;
    }


private:
    static void rhs( const state_type &x , state_type &dxdt , const double t )
    { 
	dxdt[0] = pow( t, p );
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
    std::cout << boost::format( "%-20s%-20s%-20s\n" )
	% "Estimated order" % "defined order" % "first significant error";
}

BOOST_AUTO_TEST_CASE_TEMPLATE( runge_kutta_test , Stepper, runge_kutta_steppers )
{
    integrate_const_test< Stepper > tester;
    tester();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( adams_moultion_test , Stepper, adams_steppers )
{
    integrate_const_test_initialize< Stepper > tester;
    tester();
}

BOOST_AUTO_TEST_SUITE_END()
