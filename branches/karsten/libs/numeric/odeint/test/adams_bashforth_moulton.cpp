/* Boost check_implicit_euler.cpp test file

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 This file tests the use of the euler stepper

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#define BOOST_TEST_MODULE odeint_adams_bashforth_moulton

#include <utility>

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/stepper/detail/adams_bashforth_coefficients.hpp>
#include <boost/numeric/odeint/stepper/adams_bashforth.hpp>
#include <boost/numeric/odeint/stepper/adams_bashforth_moulton.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

typedef double value_type;

BOOST_AUTO_TEST_SUITE( adams_bashforth_moulton_test )

BOOST_AUTO_TEST_CASE( test_adams_bashforth_coefficients )
{
	detail::adams_bashforth_coefficients< value_type , 1 > c1;
	detail::adams_bashforth_coefficients< value_type , 2 > c2;
	detail::adams_bashforth_coefficients< value_type , 3 > c3;
	detail::adams_bashforth_coefficients< value_type , 4 > c4;
	detail::adams_bashforth_coefficients< value_type , 5 > c5;
	detail::adams_bashforth_coefficients< value_type , 6 > c6;
	detail::adams_bashforth_coefficients< value_type , 7 > c7;
	detail::adams_bashforth_coefficients< value_type , 8 > c8;
}


BOOST_AUTO_TEST_SUITE_END()
