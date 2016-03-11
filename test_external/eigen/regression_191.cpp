/*
 [auto_generated]
 integrate.cpp

 [begin_description]
 tba.
 [end_description]

 Copyright 2016 Karsten Ahnert
 Copyright 2016 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#define BOOST_TEST_MODULE odeint_array_of_eigen

#include <boost/numeric/odeint.hpp>
#include <Eigen/Core>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE( array_of_eigen )

BOOST_AUTO_TEST_CASE( array_of_eigen_case )
{
    using state_t = std::array<Eigen::Vector3d, 1>;
    //boost::numeric::odeint::runge_kutta4<state_t> stepper;
    boost::numeric::odeint::bulirsch_stoer<state_t> stepper;

    auto tracer = [](const state_t& x, state_t& dxdt, const double /*t*/){
        dxdt[0][0] = dxdt[0][1] = dxdt[0][2] = 1;
    };

    state_t position{{{2, 2, 0}}};
    double time = 0, time_step = 0.01;

    //stepper.do_step(tracer, position, time, time_step);
    stepper.try_step(tracer, position, time, time_step);
}

BOOST_AUTO_TEST_SUITE_END()