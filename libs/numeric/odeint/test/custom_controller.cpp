/*
 [auto_generated]
 libs/numeric/odeint/test/custom_controller.cpp

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

#define BOOST_TEST_MODULE odeint_custom_controller

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/stepper/custom_controller.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>

#include "dummy_odes.hpp"

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

struct find_zero_x
{
    template< class StateIn , class StateOut , class Time >
    custom_controller_result operator()( const StateIn &in , const StateOut &out , Time t , Time dt ) const
    {
        if( out[0] >= 0.0 ) return cc_success;
        else return cc_stop;
    }
};

struct find_zero_x_precise
{
    template< class StateIn , class StateOut , class Time >
    custom_controller_result operator()( const StateIn &in , const StateOut &out , Time t , Time dt ) const
    {
        if( out[0] >= 0.0 ) return cc_success;
        else
        {
            if( out[0] > -1.0e-7 ) return cc_stop;
            else return cc_fail;
        }
    }
};


struct half_step_size
{
    template< class StateIn , class StateOut , class Time >
    Time operator()( const StateIn &in , const StateOut &out , Time t , Time &dt ) const
    {
        return dt * 0.5;
        
    }
};

struct write_to_cout
{
    template< class State , class Time >
    void operator()( const State &x , Time t ) const
    {
        std::cout << t << "\t" << x[0] << "\t" << x[1] << "\n";
    }
};


BOOST_AUTO_TEST_SUITE( custom_controller_test )

BOOST_AUTO_TEST_CASE( test_rk4 )
{
    typedef boost::array< double , 2 > state_type ;
    state_type x = {{ 1.0 , 0.0 }};
    
    custom_controller< find_zero_x_precise , half_step_size , runge_kutta4< state_type > > controller;

    size_t count = 0;
    double t = 0.0;
    double dt = 0.01;
    custom_controller_result res = cc_success;
    while( ( res != cc_stop ) && ( count < 1000 ) )
    {
        std::cout << count << "\t" << t << "\t" << dt << "\t" << x[0] << "\t" << x[1] << "\n";
        res = controller.try_step( harmonic_oscillator() , x , t , dt );
        ++count;
    }

    // integrate_cond( runge_kutta4< state_type >() , harmonic_oscillator() , x , 0.0 , 100.0 , 0.01 , write_to_cout() ,
    //                 find_zero_x_precise() , half_step_size() );

    // integrate_adaptive( make_custom_controller( find_zero_x() , half_step_size() , runge_kutta4< state_type >() ) ,
    //                     harmonic_oscillator() , x , 0.0 , 100.0 , 0.01 , write_to_cout() );
}

BOOST_AUTO_TEST_SUITE_END()
