/*
 * integrate.cpp
 *
 *  Created on: Aug 17, 2011
 *      Author: mario
 */

#define BOOST_TEST_MODULE odeint_integrate_functions

#include <vector>
#include <cmath>
#include <iostream>

#include <boost/numeric/odeint/config.hpp>

#include <boost/array.hpp>
#include <boost/ref.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <boost/test/unit_test.hpp>

#include <boost/mpl/vector.hpp>

#include <boost/numeric/odeint.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;
namespace mpl = boost::mpl;


typedef double value_type;
typedef std::vector< value_type > state_type;

void lorenz( const state_type &x , state_type &dxdt , const value_type t )
{
    const value_type sigma( 10.0 );
    const value_type R( 28.0 );
    const value_type b( value_type( 8.0 ) / value_type( 3.0 ) );

    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = -b * x[2] + x[0] * x[1];
}

struct push_back_time
{
    std::vector< double >& m_times;

    push_back_time( std::vector< double > &times )
    :  m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_times.push_back( t );
    }
};

template< class Stepper >
struct perform_integrate_const_test
{
    void operator()( void )
    {
        state_type x( 3 , 10.0 );
        const value_type dt = 0.03;
        const value_type t_end = 10.0;

        std::vector< value_type > times;

        integrate_const( Stepper() , lorenz , x , 0.0 , t_end ,
                                        dt , push_back_time( times ) );

        BOOST_CHECK_EQUAL( static_cast<int>(times.size()) , static_cast<int>(floor(t_end/dt))+1 );

        for( size_t i=0 ; i<times.size() ; ++i )
        {
            //std::cout << i << std::endl;
            // check if observer was called at times 0,1,2,...
            BOOST_CHECK_SMALL( times[i] - static_cast< value_type >(i)*dt , (i+1) * 2E-16 );
        }
    }
};

template< class Stepper >
struct perform_integrate_adaptive_test
{
    void operator()( void )
    {
        state_type x( 3 , 10.0 );
        const value_type dt = 0.03;
        const value_type t_end = 10.0;

        std::vector< value_type > times;

        size_t steps = integrate_adaptive( Stepper() , lorenz , x , 0.0 , t_end ,
                                        dt , push_back_time( times ) );

        BOOST_CHECK_EQUAL( times.size() , steps+1 );

        BOOST_CHECK_SMALL( times[0] - 0.0 , 2E-16 );
        BOOST_CHECK_SMALL( times[times.size()-1] - t_end , times.size() * 2E-16 );
    }
};


template< class Stepper >
struct perform_integrate_times_test
{
    void operator()( void )
    {
        state_type x( 3 );
        x[0] = x[1] = x[2] = 10.0;

        const value_type dt = 0.03;

        std::vector< double > times;

        // simple stepper
        integrate_times( Stepper() , lorenz , x , boost::counting_iterator<int>(0) , boost::counting_iterator<int>(10) ,
                    dt , push_back_time( times ) );

        BOOST_CHECK_EQUAL( static_cast<int>(times.size()) , 10 );

        for( size_t i=0 ; i<times.size() ; ++i )
            // check if observer was called at times 0,1,2,...
            BOOST_CHECK_EQUAL( times[i] , static_cast<double>(i) );
    }
};

template< class Stepper >
struct perform_integrate_n_steps_test
{
    void operator()( void )
    {
        state_type x( 3 );
        x[0] = x[1] = x[2] = 10.0;

        const value_type dt = 0.03;
        const int n = 200;

        std::vector< double > times;

        // simple stepper
        value_type end_time = integrate_n_steps( Stepper() , lorenz , x , 0.0 , dt , n , push_back_time( times ) );

        BOOST_CHECK_SMALL( end_time - n*dt , 2E-16 );
        BOOST_CHECK_EQUAL( static_cast<int>(times.size()) , n+1 );

        for( size_t i=0 ; i<times.size() ; ++i )
            // check if observer was called at times 0,1,2,...
            BOOST_CHECK_SMALL( times[i] - static_cast< value_type >(i)*dt , (i+1) * 2E-16 );
    }
};



class stepper_methods : public mpl::vector<
    euler< state_type > ,
    modified_midpoint< state_type > ,
    runge_kutta4< state_type > ,
    runge_kutta_cash_karp54< state_type > ,
    runge_kutta_dopri5< state_type > ,
    runge_kutta_fehlberg78< state_type > ,
    controlled_runge_kutta< runge_kutta_cash_karp54< state_type > > ,
    controlled_runge_kutta< runge_kutta_dopri5< state_type > > ,
    controlled_runge_kutta< runge_kutta_fehlberg78< state_type > > ,
    bulirsch_stoer< state_type > ,
    dense_output_runge_kutta< controlled_runge_kutta< runge_kutta_dopri5< state_type > > > ,
    bulirsch_stoer_dense_out< state_type >
> { };



BOOST_AUTO_TEST_SUITE( integrate_test )

BOOST_AUTO_TEST_CASE_TEMPLATE( integrate_const_test_case , Stepper, stepper_methods )
{
    perform_integrate_const_test< Stepper > tester;
    tester();
}


BOOST_AUTO_TEST_CASE_TEMPLATE( integrate_adaptive_test_case , Stepper, stepper_methods )
{
    perform_integrate_adaptive_test< Stepper > tester;
    tester();
}


BOOST_AUTO_TEST_CASE_TEMPLATE( integrate_times_test_case , Stepper, stepper_methods )
{
    perform_integrate_times_test< Stepper > tester;
    tester();
}

class simple_stepper_methods : public mpl::vector<
    euler< state_type > ,
    modified_midpoint< state_type > ,
    runge_kutta4< state_type > ,
    runge_kutta_cash_karp54< state_type > ,
    runge_kutta_dopri5< state_type > ,
    runge_kutta_fehlberg78< state_type >
> { };

BOOST_AUTO_TEST_CASE_TEMPLATE( integrate_n_steps_test_case , Stepper, simple_stepper_methods )
{
    perform_integrate_n_steps_test< Stepper > tester;
    tester();
}

BOOST_AUTO_TEST_SUITE_END()
