/* Boost check_implicit_euler.cpp test file

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 This file tests the use of the euler stepper

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#define BOOST_TEST_MODULE odeint_implicit_euler

#include <utility>

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/stepper/implicit_euler.hpp>
#include <boost/numeric/odeint/algebra/external/ublas_resize.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

typedef double value_type;
typedef boost::numeric::ublas::vector< value_type > state_type;
typedef boost::numeric::ublas::matrix< value_type > matrix_type;


void sys( const state_type &x , state_type &dxdt , const value_type t )
{
    dxdt( 0 ) = x( 0 ) + 2 * x( 1 );
    dxdt( 1 ) = x( 1 );
}

void jacobi( const state_type &x , matrix_type &jacobi , const value_type t )
{
    jacobi( 0 , 0 ) = 1;
    jacobi( 0 , 1 ) = 2;
    jacobi( 1 , 0 ) = 0;
    jacobi( 1 , 1 ) = 1;
}

BOOST_AUTO_TEST_SUITE( implicit_euler_test )

BOOST_AUTO_TEST_CASE( test_euler )
{
    implicit_euler< value_type > stepper;
    state_type x( 2 );
    x(0) = 0.0; x(1) = 1.0;

    value_type eps = 1E-12;

    stepper.do_step( std::make_pair( sys , jacobi ) , x , 0.0 , 0.1 );

    using std::abs;

    // compare with analytic solution of above system
    BOOST_CHECK_MESSAGE( abs( x(0) - 20.0/81.0 ) < eps , x(0) - 20.0/81.0 );
    BOOST_CHECK_MESSAGE( abs( x(1) - 10.0/9.0 ) < eps , x(0) - 10.0/9.0 );

}

BOOST_AUTO_TEST_SUITE_END()
