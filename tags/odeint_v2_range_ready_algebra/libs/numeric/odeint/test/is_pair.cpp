/*
 * is_pair.cpp
 *
 *  Created on: Feb 12, 2011
 *      Author: karsten
 */

#define BOOST_TEST_MODULE odeint_standard_algebra

#include <utility>

#include <boost/test/unit_test.hpp>
#include <boost/static_assert.hpp>

#include <boost/numeric/odeint/stepper/detail/is_pair.hpp>

using namespace boost::numeric::odeint;



BOOST_AUTO_TEST_SUITE( is_pair_test )

BOOST_AUTO_TEST_CASE( is_pair )
{
	typedef std::pair< int , int > type1;
	typedef std::pair< int& , int > type2;
	typedef std::pair< int , int& > type3;
	typedef std::pair< int& , int& > type4;
	typedef std::pair< const int , int > type5;
	typedef std::pair< const int& , int > type6;

	BOOST_STATIC_ASSERT(( detail::is_pair< type1 >::value ));
	BOOST_STATIC_ASSERT(( detail::is_pair< type2 >::value ));
	BOOST_STATIC_ASSERT(( detail::is_pair< type3 >::value ));
	BOOST_STATIC_ASSERT(( detail::is_pair< type4 >::value ));
	BOOST_STATIC_ASSERT(( detail::is_pair< type5 >::value ));
	BOOST_STATIC_ASSERT(( detail::is_pair< type6 >::value ));
}

BOOST_AUTO_TEST_SUITE_END()
