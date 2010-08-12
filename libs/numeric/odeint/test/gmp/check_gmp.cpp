/* Boost check_gmp.cpp test file

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 This file tests the odeint library with the gmp arbitrary precision types

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#include <vector>

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;


//#define BOOST_TEST_DYN_LINK

void test_gmp()
{
    //BOOST_CHECK( 1 == 0 );
}

test_suite* init_unit_test_suite( int argc , char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE( "check gmp" );

    test->add( BOOST_TEST_CASE(test_gmp) );

    return test;
}

