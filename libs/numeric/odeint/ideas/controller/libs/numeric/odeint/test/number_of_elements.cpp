/*
 * number_of_elements.cpp
 *
 *  Created on: Feb 12, 2011
 *      Author: karsten
 */

#define BOOST_TEST_MODULE odeint_number_of_elements

#include <vector>
#include <list>
#include <boost/array.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/util/number_of_elements.hpp>



using namespace boost::numeric::odeint;



BOOST_AUTO_TEST_SUITE( number_of_elements_test )

BOOST_AUTO_TEST_CASE( test_number_of_elements_vector )
{
    typedef std::vector< double > state_type;
    state_type v( 10 );
    BOOST_CHECK_EQUAL( size_t( 10 ) , number_of_elements( v ) );
}

BOOST_AUTO_TEST_CASE( test_number_of_elements_boost_array )
{
    typedef boost::array< double , 5 > state_type;
    state_type v;
    BOOST_CHECK_EQUAL( size_t( 5 ) , number_of_elements( v ) );

}

BOOST_AUTO_TEST_CASE( test_number_of_elements_list )
{
    typedef std::list< double > state_type;
    state_type v( 10 );
    BOOST_CHECK_EQUAL( size_t( 10 ) , number_of_elements( v ) );
}

BOOST_AUTO_TEST_CASE( test_number_of_elements_ublas_vector )
{
    typedef boost::numeric::ublas::vector< double > state_type;
    state_type v( 10 );
    BOOST_CHECK_EQUAL( size_t( 10 ) , number_of_elements( v ) );

}

BOOST_AUTO_TEST_SUITE_END()
