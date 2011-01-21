/*
 * check_operations.cpp
 *
 *  Created on: Jan 20, 2011
 *      Author: karsten
 */

#include <cmath>
#include <complex>
#include <utility>
#include <tr1/array>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <boost/range.hpp>

#include <boost/numeric/odeint/algebra/standard_operations.hpp>
#include <boost/numeric/odeint/algebra/standard_algebra.hpp>

using boost::numeric::odeint::standard_operations;
using boost::numeric::odeint::standard_algebra;



BOOST_AUTO_TEST_SUITE( standard_algebra_test )

BOOST_AUTO_TEST_CASE( test_for_each2 )
{
	std::tr1::array< double , 2 > x1 = {{ 1.0 , 1.0 }} , x2 = {{ 2.0 , 2.0 }};
	standard_algebra::for_each2( x1 , x2 , standard_operations::scale_sum1<>( 1.0 ) );
	BOOST_CHECK_CLOSE( x1[0] , 2.0 , 1.0e-10 );
	BOOST_CHECK_CLOSE( x1[1] , 2.0 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( test_for_each3 )
{
	std::tr1::array< double , 2 > x1 = {{ 1.0 , 1.0 }} , x2 = {{ 2.0 , 2.0 }} , x3 = {{ 3.0 , 3.0 }};
	standard_algebra::for_each3( x1 , x2 , x3 , standard_operations::scale_sum2<>( 1.0 , 2.0 ) );
	BOOST_CHECK_CLOSE( x1[0] , 8.0 , 1.0e-10 );
	BOOST_CHECK_CLOSE( x1[1] , 8.0 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( test_for_each4 )
{
	std::tr1::array< double , 2 > x1 = {{ 1.0 , 1.0 }} , x2 = {{ 2.0 , 2.0 }} , x3 = {{ 3.0 , 3.0 }} , x4 = {{ 4.0 , 4.0 }};
	standard_algebra::for_each4( x1 , x2 , x3 , x4 , standard_operations::scale_sum3<>( 1.0 , 2.0 , 3.0 ) );
	BOOST_CHECK_CLOSE( x1[0] , 20.0 , 1.0e-10 );
	BOOST_CHECK_CLOSE( x1[1] , 20.0 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( test_for_each5 )
{
	std::tr1::array< double , 2 > x1 = {{ 1.0 , 1.0 }} , x2 = {{ 2.0 , 2.0 }} , x3 = {{ 3.0 , 3.0 }} , x4 = {{ 4.0 , 4.0 }} , x5 = {{ 5.0 , 5.0 }};
	standard_algebra::for_each5( x1 , x2 , x3 , x4 , x5 , standard_operations::scale_sum4<>( 1.0 , 2.0 , 3.0 , 4.0 ) );
	BOOST_CHECK_CLOSE( x1[0] , 40.0 , 1.0e-10 );
	BOOST_CHECK_CLOSE( x1[1] , 40.0 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( test_for_each6 )
{
	std::tr1::array< double , 2 > x1 = {{ 1.0 , 1.0 }} , x2 = {{ 2.0 , 2.0 }} , x3 = {{ 3.0 , 3.0 }} , x4 = {{ 4.0 , 4.0 }} , x5 = {{ 5.0 , 5.0 }} , x6 = {{ 6.0 , 6.0 }};
	standard_algebra::for_each6( x1 , x2 , x3 , x4 , x5 , x6 ,standard_operations::scale_sum5<>( 1.0 , 2.0 , 3.0 , 4.0 , 5.0 ) );
	BOOST_CHECK_CLOSE( x1[0] , 70.0 , 1.0e-10 );
	BOOST_CHECK_CLOSE( x1[1] , 70.0 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( test_for_each7 )
{
	std::tr1::array< double , 2 > x1 = {{ 1.0 , 1.0 }} , x2 = {{ 2.0 , 2.0 }} , x3 = {{ 3.0 , 3.0 }} , x4 = {{ 4.0 , 4.0 }} , x5 = {{ 5.0 , 5.0 }} , x6 = {{ 6.0 , 6.0 }} , x7 = {{ 7.0 , 7.0 }};
	standard_algebra::for_each7( x1 , x2 , x3 , x4 , x5 , x6 , x7 , standard_operations::scale_sum6<>( 1.0 , 2.0 , 3.0 , 4.0 , 5.0 , 6.0 ) );
	BOOST_CHECK_CLOSE( x1[0] , 112.0 , 1.0e-10 );
	BOOST_CHECK_CLOSE( x1[1] , 112.0 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( test_for_each8 )
{
	std::tr1::array< double , 2 > x1 = {{ 1.0 , 1.0 }} , x2 = {{ 2.0 , 2.0 }} , x3 = {{ 3.0 , 3.0 }} , x4 = {{ 4.0 , 4.0 }} , x5 = {{ 5.0 , 5.0 }} , x6 = {{ 6.0 , 6.0 }} , x7 = {{ 7.0 , 7.0 }} , x8 = {{ 8.0 , 8.0 }};
	standard_algebra::for_each8( x1 , x2 , x3 , x4 , x5 , x6 , x7 , x8 , standard_operations::scale_sum7<>( 1.0 , 2.0 , 3.0 , 4.0 , 5.0 , 6.0 , 7.0 ) );
	BOOST_CHECK_CLOSE( x1[0] , 168.0 , 1.0e-10 );
	BOOST_CHECK_CLOSE( x1[1] , 168.0 , 1.0e-10 );
}


BOOST_AUTO_TEST_CASE( test_for_each2_with_range )
{
//	std::tr1::array< double , 2 > x1 = {{ 1.0 , 1.0 }};
//	std::tr1::array< double , 4 > x2 = {{ 2.0 , 3.0 , 4.0 , 5.0 }};
//	standard_algebra::for_each2( x1 , std::make_pair( x2.begin() + 1 , x2.begin() + 3 ) , standard_operations::scale_sum1<>( 1.0 ) );
//	BOOST_CHECK_CLOSE( x1[0] , 3.0 , 1.0e-10 );
//	BOOST_CHECK_CLOSE( x1[1] , 4.0 , 1.0e-10 );
}






BOOST_AUTO_TEST_SUITE_END()
