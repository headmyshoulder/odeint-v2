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



using namespace boost::unit_test;
using namespace boost::numeric::odeint;


BOOST_AUTO_TEST_SUITE( velocity_verlet_test )

BOOST_AUTO_TEST_CASE( test1 )
{
    typedef boost::array< double , 2 > state_type;
    velocity_verlet< state_type > stepper;
}

BOOST_AUTO_TEST_SUITE_END()
