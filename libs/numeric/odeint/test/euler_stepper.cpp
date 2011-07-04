/*
 * euler_stepper.cpp
 *
 *  Created on: Jul 4, 2011
 *      Author: mario
 */

#define BOOST_TEST_MODULE odeint_explicit_euler

#include <boost/test/unit_test.hpp>

#include <utility>
#include <iostream>
#include <vector>

#include <boost/numeric/odeint/stepper/explicit_euler.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

typedef double value_type;
typedef std::vector< value_type > state_type;

/* use functors, because functions don't work with msvc 10, I guess this is a bug */
struct sys
{
    void operator()( const state_type &x , state_type &dxdt , const value_type t ) const
    {
        dxdt[0] = x[0] + 2 * x[1];
        dxdt[1] = x[1];
    }
};


BOOST_AUTO_TEST_SUITE( explicit_euler_test )

BOOST_AUTO_TEST_CASE( test_euler )
{
    explicit_euler< state_type > stepper;
    state_type x( 2 );
    x[0] = 0.0; x[1] = 1.0;

    const value_type eps = 1E-12;
    const value_type dt = 0.1;

    stepper.do_step( sys() , x , 0.0 , dt );

    using std::abs;

    // compare with analytic solution of above system
    BOOST_CHECK_MESSAGE( abs( x[0] - 2.0*1.0*dt ) < eps , x[0] - 2.0*1.0*dt );
    BOOST_CHECK_MESSAGE( abs( x[1] - (1.0 + dt) ) < eps , x[1] - (1.0+dt) );

}

BOOST_AUTO_TEST_SUITE_END()
