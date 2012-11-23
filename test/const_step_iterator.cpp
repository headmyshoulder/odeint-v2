/*
 [auto_generated]
 libs/numeric/odeint/test/const_step_iterator.cpp

 [begin_description]
 This file tests the const step iterator.
 [end_description]

 Copyright 2009-2012 Karsten Ahnert
 Copyright 2009-2012 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#define BOOST_TEST_MODULE odeint_const_step_iterator

#include <iterator>
#include <algorithm>
#include <vector>

#include <boost/numeric/odeint/config.hpp>
#include <boost/array.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/mpl/vector.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <boost/numeric/odeint/iterator/const_step_iterator.hpp>
#include "dummy_steppers.hpp"

namespace mpl = boost::mpl;
using namespace boost::numeric::odeint;

struct dummy_system { };

typedef dummy_stepper::state_type state_type;
typedef dummy_stepper::value_type value_type;

BOOST_AUTO_TEST_SUITE( const_step_iterator_test )

typedef mpl::vector<
    dummy_stepper
    , dummy_dense_output_stepper
    > dummy_steppers;


BOOST_AUTO_TEST_CASE_TEMPLATE( transitivity1 , Stepper , dummy_steppers )
{
    typedef const_step_iterator< Stepper , dummy_system > stepper_iterator;

    state_type x = {{ 1.0 }};
    stepper_iterator first1( Stepper() , dummy_system() , x , 2.5 , 2.0 , 0.1 );
    stepper_iterator last1( Stepper() , dummy_system() , x );
    stepper_iterator last2( Stepper() , dummy_system() , x );

    BOOST_CHECK( first1 == last1 );
    BOOST_CHECK( first1 == last2 );
    BOOST_CHECK( last1 == last2 );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( copy_algorithm , Stepper , dummy_steppers )
{
    typedef const_step_iterator< Stepper , dummy_system > stepper_iterator;
    state_type x = {{ 1.0 }};
    std::vector< state_type > res;
    stepper_iterator first( Stepper() , dummy_system() , x , 0.0 , 0.35 , 0.1 );
    stepper_iterator last( Stepper() , dummy_system() , x );

    std::copy( first , last , std::back_insert_iterator< std::vector< state_type > >( res ) );

    BOOST_CHECK_EQUAL( res.size() , size_t( 4 ) );
    BOOST_CHECK_CLOSE( res[0][0] , 1.0 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[1][0] , 1.25 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[2][0] , 1.5 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[3][0] , 1.75 , 1.0e-14 );

    BOOST_CHECK_CLOSE( x[0] , 2.0 , 1.0e-14 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( copy_algorithm_negative_time_step , Stepper , dummy_steppers )
{
    typedef const_step_iterator< Stepper , dummy_system > stepper_iterator;
    state_type x = {{ 1.0 }};
    std::vector< state_type > res;
    stepper_iterator first( Stepper() , dummy_system() , x , 0.3 , -0.05 , -0.1 );
    stepper_iterator last( Stepper() , dummy_system() , x );

    std::copy( first , last , std::back_insert_iterator< std::vector< state_type > >( res ) );

    BOOST_CHECK_EQUAL( res.size() , size_t( 4 ) );
    BOOST_CHECK_CLOSE( res[0][0] , 1.0 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[1][0] , 1.25 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[2][0] , 1.5 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[3][0] , 1.75 , 1.0e-14 );

    BOOST_CHECK_CLOSE( x[0] , 2.0 , 1.0e-14 );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( copy_algorithm_with_factory , Stepper , dummy_steppers )
{
    state_type x = {{ 1.0 }};
    std::vector< state_type > res;
    std::copy( make_const_step_iterator_begin( Stepper() , dummy_system() , x , 0.0 , 0.35 , 0.1 ) ,
               make_const_step_iterator_end( Stepper() , dummy_system() , x ) ,
               std::back_insert_iterator< std::vector< state_type > >( res ) );

    BOOST_CHECK_EQUAL( res.size() , size_t( 4 ) );
    BOOST_CHECK_CLOSE( res[0][0] , 1.0 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[1][0] , 1.25 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[2][0] , 1.5 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[3][0] , 1.75 , 1.0e-14 );

    BOOST_CHECK_CLOSE( x[0] , 2.0 , 1.0e-14 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( copy_algorithm_with_range_factory , Stepper , dummy_steppers )
{
    state_type x = {{ 1.0 }};
    std::vector< state_type > res;
    boost::range::copy( make_const_step_range( Stepper() , dummy_system() , x , 0.0 , 0.35 , 0.1 ) ,
                        std::back_insert_iterator< std::vector< state_type > >( res ) );

    BOOST_CHECK_EQUAL( res.size() , size_t( 4 ) );
    BOOST_CHECK_CLOSE( res[0][0] , 1.0 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[1][0] , 1.25 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[2][0] , 1.5 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[3][0] , 1.75 , 1.0e-14 );

    BOOST_CHECK_CLOSE( x[0] , 2.0 , 1.0e-14 );
}




BOOST_AUTO_TEST_SUITE_END()
