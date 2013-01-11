/*
 [auto_generated]
 libs/numeric/odeint/test/ublas_vector_space.cpp

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

#define BOOST_TEST_MODULE odeint_dummy

#include <boost/test/unit_test.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/ublas_adapter.hpp>

#include <iostream>

#include "dummy_odes.hpp"

using namespace boost::unit_test;
using namespace boost::numeric::odeint;




BOOST_AUTO_TEST_SUITE( ublas_vector_space_algebra_test )

BOOST_AUTO_TEST_CASE( test_stepper )
{
    typedef boost::numeric::ublas::vector< double > state_type;
    typedef runge_kutta4< state_type , double , state_type , double , vector_space_algebra > stepper_type;
    
    state_type x( 1 );
    integrate_const( stepper_type() , constant_system_functor_standard() , x , 0.0 , 10.0 , 0.01 );
}


BOOST_AUTO_TEST_CASE( test_controlled_stepper )
{
    typedef boost::numeric::ublas::vector< double > state_type;
    typedef runge_kutta_cash_karp54< state_type , double , state_type , double , vector_space_algebra > stepper_type;
    auto controlled_stepper = make_controlled( 1.0e-6 , 1.0e-6 , stepper_type() );
    
    state_type x( 1 );
    integrate_const( controlled_stepper , constant_system_functor_standard() , x , 0.0 , 10.0 , 0.01 );
}

BOOST_AUTO_TEST_CASE( test_ublas_expressions )
{
    using namespace std;
    typedef boost::numeric::ublas::vector< double > state_type;
    state_type x( 2 ) , y( 2 );
    x[ 0 ] = -1.0;
    x[ 1 ] = 2.0;
    y[ 0 ] = 5.0;
    y[ 1 ] = -2.0;
    state_type z = abs( x ) + y;
    state_type w = min( abs( x ) , y ) + x + abs( y );
    cout << z[0] << " " << z[1] << "\n";
    cout << w[0] << " " << w[1] << "\n";
}

BOOST_AUTO_TEST_SUITE_END()
