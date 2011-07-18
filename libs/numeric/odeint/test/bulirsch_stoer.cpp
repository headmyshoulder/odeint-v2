/*
 * bulirsch_stoer.cpp
 *
 *  Created on: Jul 17, 2011
 *      Author: mario
 */

// disable checked iterator warning for msvc
#include <boost/config.hpp>
#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif

#define BOOST_TEST_MODULE odeint_bulirsch_stoer

#include <utility>
#include <iostream>

#include <boost/array.hpp>

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

typedef double value_type;
typedef boost::array< value_type , 3 > state_type;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

struct lorenz
{
    template< class State , class Deriv >
    void operator()( const State &x , Deriv &dxdt , double t ) const
    {
        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = -b * x[2] + x[0] * x[1];
    }
};

BOOST_AUTO_TEST_SUITE( bulirsch_stoer_test )

BOOST_AUTO_TEST_CASE( test_bulirsch_stoer )
{
    typedef bulirsch_stoer< state_type > stepper_type;
    stepper_type stepper( 1E-9 , 1E-9 , 1.0 , 0.0 );

    state_type x;
    x[0] = 10.0 ; x[1] = 10.0 ; x[2] = 5.0;

    double dt = 0.1;

    //stepper.try_step( lorenz() , x , t , dt );

    size_t steps = integrate_adaptive( stepper , lorenz() , x , 0.0 , 10.0 , dt );

    std::cout << "required steps: " << steps << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
