/*
 [auto_generated]
 const_step_range.cpp

 [begin_description]
 tba.
 [end_description]

 Copyright 2009-2012 Karsten Ahnert
 Copyright 2009-2012 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/range/const_step_range.hpp>
#include "dummy_odes.hpp"

#include <boost/range/algorithm/for_each.hpp>

#include <boost/config.hpp>
#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif

#define BOOST_TEST_MODULE odeint_const_step_range_test

#include <boost/test/unit_test.hpp>

#include <iostream>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

struct observer
{
    template< typename State >
    void operator()( State const& x ) const
    {
        std::cout << x << std::endl;
    }
};


BOOST_AUTO_TEST_SUITE( const_step_range_test_test )

BOOST_AUTO_TEST_CASE( test_case1 )
{
    typedef boost::array< double , 1 > state_type;
    typedef runge_kutta4< state_type > stepper_type;
    typedef constant_system_functor_standard system_type;
    typedef double time_type;
    typedef const_step_range< stepper_type , state_type , system_type , time_type > range_type;
    
    stepper_type stepper;
    state_type state;
    system_type system;
    range_type r( stepper , state , system , 0.0 , 10.0 , 0.01 );
    
    boost::for_each( r , observer() );
    
    
    
    BOOST_CHECK_EQUAL( 1 , 1 );
}

BOOST_AUTO_TEST_SUITE_END()
