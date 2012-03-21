/*
 * is_pair.cpp
 *
 *  Created on: Feb 12, 2011
 *      Author: karsten
 */

#define BOOST_TEST_MODULE odeint_ref_or_value

#include <utility>

#include <boost/test/unit_test.hpp>
#include <boost/static_assert.hpp>

#include <boost/numeric/odeint/util/ref_or_value_holder.hpp>



using namespace boost::numeric::odeint;



BOOST_AUTO_TEST_SUITE( ref_or_value_holder_test )

BOOST_AUTO_TEST_CASE( test_value_holder )
{
    int a = 1;
    ref_or_value_holder< int , false > value( a );
    value.get() = 2;

    BOOST_CHECK_EQUAL( 1 , a  );
    BOOST_CHECK_EQUAL( 2 , value.get() );
}

BOOST_AUTO_TEST_CASE( test_ref_holder )
{
    int a = 1;
    ref_or_value_holder< int , true > value( a );
    value.get() = 2;

    BOOST_CHECK_EQUAL( 2 , a  );
    BOOST_CHECK_EQUAL( 2 , value.get() );
}



BOOST_AUTO_TEST_CASE( test_const_value_holder )
{
    int a = 1;
    ref_or_value_holder< const int , false > value( a );

    BOOST_CHECK_EQUAL( 1 , a  );
    BOOST_CHECK_EQUAL( 1 , value.get() );
}

BOOST_AUTO_TEST_CASE( test_const_ref_holder )
{
    int a = 1;
    ref_or_value_holder< const int , true > value( a );

    BOOST_CHECK_EQUAL( 1 , a  );
    BOOST_CHECK_EQUAL( 1 , value.get() );
}

BOOST_AUTO_TEST_SUITE_END()
