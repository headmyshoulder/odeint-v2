/*
 [auto_generated]
 libs/numeric/odeint/test/velocity_verlet.cpp

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

#define BOOST_TEST_MODULE odeint_velocity_verlet

#include <boost/numeric/odeint/stepper/velocity_verlet.hpp>

#include <boost/array.hpp>
#include <boost/test/unit_test.hpp>
#include "../../../../../../boost/boost_1_54_0/boost/concept_check.hpp"



using namespace boost::unit_test;
using namespace boost::numeric::odeint;

struct ode
{
    template< class CoorIn , class MomentumIn , class AccelerationOut >
    void operator()( const CoorIn &q , const MomentumIn &p , AccelerationOut &a ) const
    {
        a[0] = -q[0] - p[0];
        a[1] = -q[1] - p[1];
    }
};

template< class Q , class P >
void init_state( Q &q , P &p )
{
    q[0] = 1.0 ; q[1] = 0.5;
    p[0] = 2.0 ; p[1] = -1.0;
}

typedef boost::array< double , 2 > array_type;
typedef std::vector< double > vector_type;

typedef velocity_verlet< array_type > array_stepper;
typedef velocity_verlet< vector_type > vector_stepper;

BOOST_AUTO_TEST_SUITE( velocity_verlet_test )

BOOST_AUTO_TEST_CASE( test_with_array_ref )
{
    array_stepper stepper;
    array_type q , p ;
    init_state( q , p );
    stepper.do_step( ode() , std::make_pair( boost::ref( q ) , boost::ref( p ) ) , 0.0 , 0.01 );
}

BOOST_AUTO_TEST_CASE( test_with_array_pair )
{
    array_stepper stepper;
    std::pair< array_type , array_type > xxx;
    init_state( xxx.first , xxx.second );
    stepper.do_step( ode() , xxx , 0.0 , 0.01 );
}

BOOST_AUTO_TEST_CASE( test_with_vector_ref )
{
    vector_stepper stepper;
    vector_type q( 2 ) , p( 2 );
    init_state( q , p );
    stepper.do_step( ode() , std::make_pair( boost::ref( q ) , boost::ref( p ) ) , 0.0 , 0.01 );
}

BOOST_AUTO_TEST_CASE( test_with_vector_pair )
{
    vector_stepper stepper;
    std::pair< vector_type , vector_type > x;
    init_state( x.first , x.second );
    stepper.do_step( ode() , x , 0.0 , 0.01 );
}

BOOST_AUTO_TEST_CASE( test_reset )
{
}

BOOST_AUTO_TEST_CASE( test_initialize )
{
}

BOOST_AUTO_TEST_CASE( test_resize )
{
}


BOOST_AUTO_TEST_SUITE_END()
