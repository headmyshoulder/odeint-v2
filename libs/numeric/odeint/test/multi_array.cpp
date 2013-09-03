/*
 [auto_generated]
 libs/numeric/odeint/test/multi_array.cpp

 [begin_description]
 tba.
 [end_description]

 Copyright 2009-2012 Karsten Ahnert
 Copyright 2009-2012 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/config.hpp>
#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif

#define BOOST_TEST_MODULE odeint_multi_array

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/algebra/multi_array_algebra.hpp>
#include <boost/numeric/odeint/util/multi_array_adaption.hpp>


using namespace boost::unit_test;
using namespace boost::numeric::odeint;

typedef boost::multi_array< double , 1 > vector_type;
typedef boost::multi_array< double , 2 > matrix_type;
typedef boost::multi_array< double , 3 > tensor_type;


BOOST_AUTO_TEST_SUITE( multi_array_test )

BOOST_AUTO_TEST_CASE( test_resizeable )
{
    BOOST_CHECK( is_resizeable< vector_type >::value );
    BOOST_CHECK( is_resizeable< matrix_type >::value );
    BOOST_CHECK( is_resizeable< tensor_type >::value );
}

BOOST_AUTO_TEST_CASE( test_same_size_vector1 )
{
    vector_type v1( boost::extents[12] );
    vector_type v2( boost::extents[12] );
    BOOST_CHECK( same_size( v1 , v2 ) );
}

BOOST_AUTO_TEST_CASE( test_same_size_vector2 )
{
    vector_type v1( boost::extents[12] );
    vector_type v2( boost::extents[13] );
    BOOST_CHECK( !same_size( v1 , v2 ) );
}

BOOST_AUTO_TEST_CASE( test_same_size_vector3 )
{
    vector_type v1( boost::extents[12] );
    vector_type v2( boost::extents[vector_type::extent_range(-1,11)] );
    BOOST_CHECK( same_size( v1 , v2 ) );
}

BOOST_AUTO_TEST_CASE( test_same_size_vector4 )
{
    vector_type v1( boost::extents[12] );
    vector_type v2( boost::extents[vector_type::extent_range(-1,10)] );
    BOOST_CHECK( !same_size( v1 , v2 ) );
}


BOOST_AUTO_TEST_CASE( test_same_size_matrix1 )
{
    matrix_type m1( boost::extents[12][4] );
    matrix_type m2( boost::extents[12][4] );
    BOOST_CHECK( same_size( m1 , m2 ) );
}

BOOST_AUTO_TEST_CASE( test_same_size_matrix2 )
{
    matrix_type m1( boost::extents[12][4] );
    matrix_type m2( boost::extents[12][3] );
    BOOST_CHECK( !same_size( m1 , m2 ) );
}

BOOST_AUTO_TEST_CASE( test_same_size_matrix3 )
{
    matrix_type m1( boost::extents[12][matrix_type::extent_range(-1,2)] );
    matrix_type m2( boost::extents[12][3] );
    BOOST_CHECK( same_size( m1 , m2 ) );
}

BOOST_AUTO_TEST_CASE( test_same_size_matrix4 )
{
    matrix_type m1( boost::extents[12][matrix_type::extent_range(-1,1)] );
    matrix_type m2( boost::extents[12][3] );
    BOOST_CHECK( !same_size( m1 , m2 ) );
}

BOOST_AUTO_TEST_CASE( test_same_size_tensor1 )
{
    tensor_type t1( boost::extents[ tensor_type::extent_range( -2 , 9 ) ]
                                  [ tensor_type::extent_range( 5 , 11 ) ]
                                  [ tensor_type::extent_range( -1 , 2 ) ] );
    tensor_type t2( boost::extents[ tensor_type::extent_range( 2 , 13 ) ]
                                  [ tensor_type::extent_range( -1 , 5 ) ]
                                  [ tensor_type::extent_range( 1 , 4 ) ] );
    BOOST_CHECK( same_size( t1 , t2 ) );    
}

// Tests for tensor_type

BOOST_AUTO_TEST_CASE( test_resize_vector1 )
{
//     vector_type v1( boost::extents[4] );
//     vector_type v2;
//     resize( v2 , v1 );
//     BOOST_CHECK( same_size( v1 , v2 ) );
}


BOOST_AUTO_TEST_SUITE_END()
