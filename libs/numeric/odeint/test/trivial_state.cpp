/* Boost trivial_state.cpp test file

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 This file tests the use of the euler stepper

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

// disable checked iterator warning for msvc
#include <boost/config.hpp>
#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif

#define BOOST_TEST_MODULE odeint_trivial_state

#include <boost/test/unit_test.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/at.hpp>

#include <boost/numeric/odeint.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

namespace mpl = boost::mpl;

typedef double state_type;
typedef double value_type;
typedef double deriv_type;
typedef double time_type;

void constant_system( const state_type &x , state_type &dxdt , const time_type t ) { dxdt = 1.0; }


BOOST_AUTO_TEST_SUITE( trivial_state_test )

/* test different do_step methods of simple steppers */

typedef mpl::vector<
        euler< state_type , value_type , deriv_type , time_type ,
               vector_space_algebra , default_operations , never_resizer > ,
        runge_kutta4< state_type , value_type , deriv_type , time_type ,
                      vector_space_algebra , default_operations , never_resizer >
        >::type stepper_types;


BOOST_AUTO_TEST_CASE_TEMPLATE( test_do_step , T, stepper_types )
{
	typedef T stepper_type;
	stepper_type stepper;
	state_type x = 0.0;
	time_type t = 0.0;
	time_type dt = 0.1;
	stepper.do_step( constant_system , x , t , dt );

	// this overload is not allowed if the types of dxdt and dt are the same
	// deriv_type dxdt = 1.0;
	// stepper.do_step( constant_system , x , dxdt , t , dt );

	state_type x_out;
	stepper.do_step( constant_system , x , t , x_out , dt );
}


/* test integrate_adaptive with controlled steppers */

typedef mpl::vector<
    /* bigest problem: integrate_adaptive doesn't work with controlled cash_karp! */
    runge_kutta_cash_karp54< state_type , value_type , deriv_type , time_type ,
     			     vector_space_algebra , default_operations , never_resizer > ,
    /* dopri5 is working though */
    runge_kutta_dopri5< state_type , value_type , deriv_type , time_type ,
			vector_space_algebra , default_operations , never_resizer >
    > error_stepper_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( test_integrate , T , error_stepper_types )
{
    typedef T stepper_type;
    state_type x = 0.0;
    time_type t0 = 0.0;
    time_type t1 = 1.0;
    time_type dt = 0.1;
    integrate_adaptive( make_controlled< stepper_type >( 1e-6 , 1e-6 ) , constant_system , x , t0 , t1 , dt );
}

BOOST_AUTO_TEST_SUITE_END()
