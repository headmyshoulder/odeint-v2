/*
 * const_step_ode_iterator.cpp
 *
 *  Created on: Jun 26, 2012
 *      Author: karsten
 */

#define BOOST_TEST_MODULE odeint_const_step_iterator

#include <iostream>
#include <utility>
#include <algorithm>

#include <boost/test/unit_test.hpp>
#include <boost/static_assert.hpp>
#include <boost/array.hpp>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/iterator/const_step_iterator.hpp>


using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

//[ system_function_without_perturbations
struct lorenz
{
    template< class State , class Deriv >
    void operator()( const State &x_ , Deriv &dxdt_ , double t ) const
    {
        typename boost::range_iterator< const State >::type x = boost::begin( x_ );
        typename boost::range_iterator< Deriv >::type dxdt = boost::begin( dxdt_ );

        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = -b * x[2] + x[0] * x[1];
    }
};
//]




BOOST_AUTO_TEST_SUITE( const_step_iterator_test )

BOOST_AUTO_TEST_CASE( empty_test )
{
}

struct printer
{
    template< class State >
    void operator()( const State &s )
    {
        std::cerr << "huhu" << std::endl;
    }
};

BOOST_AUTO_TEST_CASE( compile_test )
{
    typedef boost::array< double , 3 > state_type;
    runge_kutta4< state_type > stepper;
    state_type x = {{ 10.0 , 10.0 , 10.0 }};
    std::for_each( make_const_step_iterator( stepper , lorenz() , x , 0.0 , 0.01 ) ,
                   make_const_step_iterator( stepper , lorenz() , x , 10.0 , 0.01 ) , printer() );
}

// BOOST_AUTO_TEST_CASE( test_is_pair )
// {
// 	typedef std::pair< int , int > type1;
// 	typedef std::pair< int& , int > type2;
// 	typedef std::pair< int , int& > type3;
// 	typedef std::pair< int& , int& > type4;
// 	typedef std::pair< const int , int > type5;
// 	typedef std::pair< const int& , int > type6;

// 	BOOST_STATIC_ASSERT(( is_pair< type1 >::value ));
// 	BOOST_STATIC_ASSERT(( is_pair< type2 >::value ));
// 	BOOST_STATIC_ASSERT(( is_pair< type3 >::value ));
// 	BOOST_STATIC_ASSERT(( is_pair< type4 >::value ));
// 	BOOST_STATIC_ASSERT(( is_pair< type5 >::value ));
// 	BOOST_STATIC_ASSERT(( is_pair< type6 >::value ));
// }

BOOST_AUTO_TEST_SUITE_END()
