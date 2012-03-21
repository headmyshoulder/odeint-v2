/*
 * is_pair.cpp
 *
 *  Created on: Feb 12, 2011
 *      Author: karsten
 */

#define BOOST_TEST_MODULE odeint_error_checker_max_norm

#include <utility>

#include <boost/test/unit_test.hpp>
#include <boost/static_assert.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/stepper/error_checker_max_norm.hpp>

#include <boost/array.hpp>

using namespace boost::numeric::odeint;



BOOST_AUTO_TEST_SUITE( error_checker_max_norm_test )

BOOST_AUTO_TEST_CASE( test_version1 )
{
    // computes err = max( x_err / ( eps_abs + eps_rel * max( x , x_old ) ) );

    range_algebra algebra;
    error_checker_max_norm< double , range_algebra , default_operations > error_checker( algebra , 1.0e-3 , 1.0e-4 );

    double dt , err;
    boost::array< double , 2 > x , x_old , x_err;

    x[0] = 1.0 ; x[1] = 2.0;
    x_old[0] = 0.5 ; x_old[1] = 1.5;
    x_err[0] = 0.01 ; x_err[1] = 0.02;

    // err1 = 0.01 / ( 0.001 + 0.0001 * max( 1.0 , 0.5 ) )
    // err1 = 0.01 / ( 0.0011 )
    // err1 = 9.090909090
    // err2 = 0.02 / ( 0.001 + 0.0001 * max( 2.0 , 1.5 ) )
    // err2 = 0.02 / ( 0.0012 )
    // err2 = 16.6666666666666666
    // err = max( err1 , err2 )


    err = error_checker.error( x_old , x , x_err , dt );
    BOOST_CHECK_CLOSE( err , 16.66666666666 , 1.0e-10 );

    x_err[1] = 0.00002;
    err = error_checker.error( x_old , x , x_err , dt );
    BOOST_CHECK_CLOSE( err , 9.0909090909091 , 1.0e-10 );

    x_err[0] *= -1.0;
    err = error_checker.error( x_old , x , x_err , dt );
    BOOST_CHECK_CLOSE( err , 9.0909090909091 , 1.0e-10 );

    x[0] *= -1.0;
    err = error_checker.error( x_old , x , x_err , dt );
    BOOST_CHECK_CLOSE( err , 9.0909090909091 , 1.0e-10 );

    x_old[0] *= -1.0;
    err = error_checker.error( x_old , x , x_err , dt );
    BOOST_CHECK_CLOSE( err , 9.0909090909091 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( test_version2 )
{
    range_algebra algebra;
    error_checker_max_norm< double , range_algebra , default_operations > error_checker( algebra , 1.0e-3 , 1.0e-4 , 0.5 , 1.5 );

    double dt = 0.02 , err;
    boost::array< double , 2 > x , x_old , x_err , dxdt_old , dxdt;

    x[0] = 1.0 ; x[1] = 2.0;
    x_old[0] = 0.5 ; x_old[1] = 1.5;
    x_err[0] = 0.01 ; x_err[1] = 0.02;
    dxdt_old[0] = 2.0 ; dxdt_old[1] = 4.0;
    dxdt[0] = 1000.0 , dxdt[1] = 255555.0;

    // err1 = 0.01 / ( 0.001 + 0.0001 * ( 0.5 * 0.5 + 0.02 * 1.5 * 2.0 ) )
    // err1 = 0.01 / ( 0.001 + 0.0001 * ( 0.25 + 0.06 ) )
    // err1 = 0.01 / ( 0.001031
    // err1 = 9.699321048
    // err2 = 0.02 / ( 0.001 + 0.0001 * ( 0.5 * 1.5 + 0.02 * 1.5 * 4.0 ) )
    // err2 = 0.02 / ( 0.001 + 0.0001 * ( 0.75 + 0.12 ) )
    // err2 = 0.02 / ( 0.001087 )
    // err2 = 18.399264029

    err = error_checker.error( x_old , x , dxdt_old , x_err , dt );
    BOOST_CHECK_CLOSE( err , 18.399264029 , 1.0e-8 );

    err = error_checker.error( x_old , x , dxdt_old , dxdt , x_err , dt );
//    BOOST_CHECK_CLOSE( err , 18.399264029 , 1.0e-8 );

    x_err[1] = 0.00002;
    err = error_checker.error( x_old , x , dxdt_old , x_err , dt );
    BOOST_CHECK_CLOSE( err , 9.699321048 , 1.0e-8 );

    x_err[0] *= -1.0;
    err = error_checker.error( x_old , x , dxdt_old , x_err , dt );
    BOOST_CHECK_CLOSE( err , 9.699321048 , 1.0e-8 );

    x_old[0] *= -1.0;
    err = error_checker.error( x_old , x , dxdt_old , x_err , dt );
    BOOST_CHECK_CLOSE( err , 9.699321048 , 1.0e-8 );

    dxdt_old[0] *= 1.0;
    err = error_checker.error( x_old , x , dxdt_old , x_err , dt );
    BOOST_CHECK_CLOSE( err , 9.699321048 , 1.0e-8 );
}



BOOST_AUTO_TEST_SUITE_END()
