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

/* generic test for all steppers that support integrate_const */
template< class Stepper >
struct integrate_const_test
{
    void operator()( int nSteps = 1 )
    {
	double tolerance = 1.0e-13;

	state_type x;
        Stepper stepper;
	const int o = stepper.order()+1; 
	for ( p = 0; p < o-1; p++ )
	    {
		Stepper stepper1;
		x[0] = 1.0;
		integrate_const( stepper1, rhs, x, 0.0, 1.0, 1.0/nSteps );
		BOOST_CHECK_LT( fabs( x[0]-1.0/(1.0+p) - 1.0 ), tolerance  );
		std::cout << fabs( x[0]-1.0/(1.0+p)-1.0 ) << std::endl;
	    }	
        std::cout << std::endl;
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
    adams_bashforth< 2, state_type >,
    adams_bashforth< 3, state_type >,
    adams_bashforth< 4, state_type >,
    adams_bashforth_moulton< 2, state_type >,
    adams_bashforth_moulton< 3, state_type >,
    adams_bashforth_moulton< 4, state_type >
    // \TODO: write tests for order bigger than 4. 
    //        initialize with fehlberg
    > adams_steppers;

typedef mpl::vector<
    symplectic_euler< state_type > ,
    symplectic_rkn_sb3a_mclachlan< state_type >
    > symplectic_steppers;


BOOST_AUTO_TEST_CASE_TEMPLATE( runge_kutta_test , Stepper, runge_kutta_steppers )
{
    integrate_const_test< Stepper > tester;
    tester();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( adams_moultion_test , Stepper, adams_steppers )
{
    integrate_const_test< Stepper > tester;
    tester( 16 );
}

// \TODO: add symplectic tests

BOOST_AUTO_TEST_SUITE_END()
