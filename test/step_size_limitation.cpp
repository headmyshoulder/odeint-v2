/*
 [auto_generated]
 libs/numeric/odeint/test/step_size_limitation.cpp

 [begin_description]
 Tests the step size limitation functionality
 [end_description]

 Copyright 2015 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#define BOOST_TEST_MODULE odeint_integrate_times

#include <boost/test/unit_test.hpp>

#include <utility>
#include <iostream>
#include <vector>

#include <boost/numeric/odeint.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

typedef double value_type;
typedef std::vector< value_type > state_type;


void damped_osc( const state_type &x , state_type &dxdt , const value_type t )
{
    const value_type gam( 0.1);

    dxdt[0] = x[1];
    dxdt[1] = -x[0] - gam*x[1];
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

BOOST_AUTO_TEST_SUITE( step_size_limitation_test )

BOOST_AUTO_TEST_CASE( test_step_adjuster )
{
    // first use adjuster without step size limitation
    default_step_adjuster<double, double> step_adjuster;
    const double dt = 0.1;
    double dt_new = step_adjuster.decrease_step(dt, 1.5, 2);
    BOOST_CHECK(dt_new < dt*2.0/3.0);

    dt_new = step_adjuster.increase_step(dt, 0.8, 1);
    // for errors > 0.5 no increase is performed
    BOOST_CHECK(dt_new == dt);

    dt_new = step_adjuster.increase_step(dt, 0.4, 1);
    // smaller errors should lead to step size increase
    std::cout << dt_new << std::endl;
    BOOST_CHECK(dt_new > dt);


    // now test with step size limitation max_dt = 0.1
    default_step_adjuster<double, double>
        limited_adjuster(dt);

    dt_new = limited_adjuster.decrease_step(dt, 1.5, 2);
    // decreasing works as before
    BOOST_CHECK(dt_new < dt*2.0/3.0);

    dt_new = limited_adjuster.decrease_step(2*dt, 1.1, 2);
    // decreasing a large step size should give max_dt
    BOOST_CHECK(dt_new == dt);

    dt_new = limited_adjuster.increase_step(dt, 0.8, 1);
    // for errors > 0.5 no increase is performed, still valid
    BOOST_CHECK(dt_new == dt);

    dt_new = limited_adjuster.increase_step(dt, 0.4, 1);
    // but even for smaller errors, we should at most get 0.1
    BOOST_CHECK(dt_new == dt);

    dt_new = limited_adjuster.increase_step(0.9*dt, 0.1, 1);
    std::cout << dt_new << std::endl;
    // check that we don't increase beyond max_dt
    BOOST_CHECK(dt_new == dt);
}


template<class ControlledStepper>
void test_controlled_stepper(ControlledStepper stepper, const double max_dt)
{
    state_type x( 2 );
    x[0] = x[1] = 10.0;

    std::vector<double> times;

    integrate_adaptive(stepper, damped_osc, x, 0.0, 100*max_dt, max_dt, push_back_time(times));
    // check that dt remains at exactly max_dt

    std::cout << "STEP: " << times.size() << " , END: " << times[times.size()-1] << std::endl;

    for( size_t i=0 ; i<times.size() ; ++i )
        // check if observer was called at times 0,1,2,...
        BOOST_CHECK_SMALL( times[i] - static_cast<double>(i)*max_dt , 1E-15);
    times.clear();

    return;
    // this should also work when we provide some bigger initial dt
    integrate_adaptive(stepper, damped_osc, x, 0.0, 100*max_dt, 10*max_dt, push_back_time(times));
    // check that dt remains at exactly max_dt
    for( size_t i=0 ; i<times.size() ; ++i )
        // check if observer was called at times 0,1,2,...
        BOOST_CHECK_SMALL( times[i] - static_cast<double>(i)*max_dt , 1E-15);
    times.clear();
}


BOOST_AUTO_TEST_CASE( test_controlled )
{
    const double max_dt = 0.01;

    test_controlled_stepper(make_controlled(1E-2, 1E-2, max_dt,
                                            runge_kutta_dopri5<state_type>()),
                            max_dt);
    test_controlled_stepper(make_controlled(1E-2, 1E-2, -max_dt,
                                            runge_kutta_dopri5<state_type>()),
                            -max_dt);
    std::cout << "DOPRI5 FINISHED" << std::endl;

    test_controlled_stepper(make_controlled(1E-2, 1E-2, max_dt,
                                            runge_kutta_cash_karp54<state_type>()),
                            max_dt);
    test_controlled_stepper(make_controlled(1E-2, 1E-2, -max_dt,
                                            runge_kutta_cash_karp54<state_type>()),
                            -max_dt);
    std::cout << "RK-CK54 FINISHED" << std::endl;

    test_controlled_stepper(bulirsch_stoer<state_type>(1E-2, 1E-2, 1.0, 1.0, max_dt),
                            max_dt);
    test_controlled_stepper(bulirsch_stoer<state_type>(1E-2, 1E-2, 1.0, 1.0, -max_dt),
                            -max_dt);
    std::cout << "BS FINISHED" << std::endl;
}



BOOST_AUTO_TEST_CASE( test_dense_out )
{
    state_type x( 2 );
    x[0] = x[1] = 10.0;
    const double max_dt = 0.01;

    std::vector<double> times;

    integrate_adaptive(make_dense_output(1E-2, 1E-2, max_dt, runge_kutta_dopri5<state_type>()),
                       damped_osc, x, 0.0, 1.0, max_dt, push_back_time(times));
    // check that dt remains at exactly max_dt
    for( size_t i=0 ; i<times.size() ; ++i )
        // check if observer was called at times 0,1,2,...
        BOOST_CHECK_SMALL( times[i] - static_cast<double>(i)*max_dt , 1E-15);
    times.clear();

    // this should also work when we provide some bigger initial dt
    integrate_adaptive(make_dense_output(1E-2, 1E-2, max_dt, runge_kutta_dopri5<state_type>()),
                       damped_osc, x, 0.0, 1.0, 0.1, push_back_time(times));
    // check that dt remains at exactly max_dt
    for( size_t i=0 ; i<times.size() ; ++i )
        // check if observer was called at times 0,1,2,...
        BOOST_CHECK_SMALL( times[i] - static_cast<double>(i)*max_dt , 1E-15);
    times.clear();
}

BOOST_AUTO_TEST_SUITE_END()
