/* Boost numeric test of the runge kutta steppers test file

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 This file tests the use of the adams bashforth stepper

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

// disable checked iterator warning for msvc
#include <boost/config.hpp>
#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif

#define BOOST_TEST_MODULE numeric_runge_kutta

#include <iostream>
#include <cmath>

#include <boost/array.hpp>

#include <boost/test/unit_test.hpp>

#include <boost/mpl/list.hpp>
#include <boost/mpl/size_t.hpp>
#include <boost/mpl/range_c.hpp>


#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

typedef double value_type;

struct osc
{
	template< class State , class Deriv , class Value >
	void operator()( const State &x , Deriv &dxdt , const Value &t ) const
	{
            dxdt[0] = x[1];
            dxdt[1] = -x[0];
	}
};


BOOST_AUTO_TEST_SUITE( numeric_runge_kutta_test )

BOOST_AUTO_TEST_CASE( test_order )
{
    typedef boost::array< double , 2 > state_type;
    
    runge_kutta4< state_type > stepper;
    const int o = stepper.order()+1; //order of the error is order of approximation + 1

    const state_type x0 = {{ 0.0 , 1.0 }};
    state_type x1;
    const double t = 0.0;
    /* do a first step with dt=0.1 to get an estimate on the prefactor of the error dx = f * dt**(order+1) */
    double dt = 0.1;
    stepper.do_step( osc() , x0 , t , x1 , dt );
    const double f = 2.0 * std::abs( sin(dt) - x1[0] ) / std::pow( dt , o );

    std::cout << f << std::endl;

    /* as long as we have errors above machine precision */
    while( std::pow( dt , o ) > 1E-16 )
    {
        dt *= 0.5;
        stepper.do_step( osc() , x0 , t , x1 , dt );
        std::cout << "Testing dt=" << dt << std::endl;
        BOOST_CHECK_SMALL( std::abs( sin(dt) - x1[0] ) , f*std::pow( dt , o ) );
    }
}

BOOST_AUTO_TEST_SUITE_END()
