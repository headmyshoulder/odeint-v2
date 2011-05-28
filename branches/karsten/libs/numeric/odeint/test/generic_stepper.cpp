/* Boost generic_stepper.cpp test file

 Copyright 2011 Karsten Ahnert
 Copyright 2011 Mario Mulansky

 This file tests the use of the generic stepper

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#define BOOST_TEST_MODULE odeint_generic_stepper

#include <utility>

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/stepper/explicit_generic_rk.hpp>
#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>

#include <boost/array.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

typedef double value_type;
typedef boost::array< value_type , 2 > state_type;

void sys( const state_type &x , state_type &dxdt , const value_type &t )
{
    dxdt[ 0 ] = x[ 0 ] + 2 * x[ 1 ];
    dxdt[ 1 ] = x[ 1 ];
}

typedef explicit_generic_rk< 4 , 4 , state_type> stepper_type;

const boost::array< double , 1 > a1 = {{ 0.5 }};
const boost::array< double , 2 > a2 = {{ 0.0 , 0.5 }};
const boost::array< double , 3 > a3 = {{ 0.0 , 0.0 , 1.0 }};

const stepper_type::coef_a_type a = fusion::make_vector( a1 , a2 , a3 );
const stepper_type::coef_b_type b = {{ 1.0/6 , 1.0/3 , 1.0/3 , 1.0/6 }};
const stepper_type::coef_c_type c = {{ 0.0 , 0.5 , 0.5 , 1.0 }};

typedef explicit_rk4< state_type > rk4_stepper_type;

BOOST_AUTO_TEST_SUITE( generic_stepper_test )

BOOST_AUTO_TEST_CASE( test_generic_stepper )
{
	stepper_type stepper( a , b , c );

	rk4_stepper_type rk4;

	typedef stepper_type::state_type state_type;
	typedef stepper_type::value_type stepper_value_type;
	typedef stepper_type::deriv_type deriv_type;
	typedef stepper_type::time_type time_type;

	state_type x = {{ 0.0 , 1.0 }};
	state_type y = x;
	
	stepper.do_step( sys , x , 0.0 , 0.1 );

	rk4.do_step( sys , y , 0.0 , 0.1 );

	// compare with analytic solution of above system
	BOOST_CHECK_EQUAL( x[0] , y[0] );
	BOOST_CHECK_EQUAL( x[1] , y[1] );

}

BOOST_AUTO_TEST_SUITE_END()