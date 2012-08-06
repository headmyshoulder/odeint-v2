/*
 * const_step_iterator.cpp
 *
 *  Created on: Aug 6, 2012
 *      Author: Karsten Ahnert
 */

#define BOOST_TEST_MODULE odeint_const_step_iterator

#include <iterator>
#include <algorithm>
#include <vector>

#include <boost/numeric/odeint/config.hpp>
#include <boost/array.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <boost/numeric/odeint/iterator/const_step_iterator.hpp>


using namespace boost::numeric::odeint;

struct stepper_mock
{
    typedef double value_type;
    typedef value_type time_type;
    typedef boost::array< value_type , 1 > state_type;
    typedef state_type deriv_type;
    typedef unsigned short order_type;
    typedef stepper_tag stepper_category;

    order_type order( void ) const { return 1; }

    template< class System >
    void do_step( System sys , state_type &x , time_type t , time_type dt ) const
    {
        x[0] += 0.25;
    }
};

struct dummy_system { };

typedef const_step_iterator< stepper_mock , dummy_system > stepper_iterator;
typedef stepper_mock::state_type state_type;
typedef stepper_mock::value_type value_type;

BOOST_AUTO_TEST_SUITE( const_step_iterator_test )

BOOST_AUTO_TEST_CASE( transitivity1 )
{
    state_type x = {{ 1.0 }};
    stepper_iterator first1( stepper_mock() , dummy_system() , x , 1.5 , 0.1 , true );
    stepper_iterator first2( stepper_mock() , dummy_system() , x , 2.0 , 0.1 , false );
    stepper_iterator last( stepper_mock() , dummy_system() , x , 1.0 , 0.1 , false );

    BOOST_CHECK( first1 == last );
    BOOST_CHECK( first2 == last );

    // this one fails
    BOOST_CHECK( first1 == first2 );
}

BOOST_AUTO_TEST_CASE( transitivity2 )
{
    state_type x = {{ 1.0 }};
    stepper_iterator first1( stepper_mock() , dummy_system() , x , 1.5 , 0.1 , true );
    stepper_iterator first2( stepper_mock() , dummy_system() , x , 2.0 , 0.1 , false );
    stepper_iterator last( stepper_mock() , dummy_system() , x , 1.0 , 0.1 , false );
}


BOOST_AUTO_TEST_CASE( copy_algo )
{
    state_type x = {{ 1.0 }};
    std::vector< state_type > res;
    stepper_iterator first( stepper_mock() , dummy_system() , x , 0.0 , 0.1 , true );
    stepper_iterator last( stepper_mock() , dummy_system() , x , 0.35 , 0.1 , false );

    std::copy( first , last , std::back_insert_iterator< std::vector< state_type > >( res ) );

    BOOST_CHECK_EQUAL( res.size() , size_t( 4 ) );
    BOOST_CHECK_CLOSE( res[0][0] , 1.0 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[1][0] , 1.25 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[2][0] , 1.5 , 1.0e-14 );
    BOOST_CHECK_CLOSE( res[3][0] , 1.75 , 1.0e-14 );
     
}



BOOST_AUTO_TEST_SUITE_END()
