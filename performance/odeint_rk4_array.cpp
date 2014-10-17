/*
 * odeint_rk4_lorenz_def_alg.cpp
 *
 * Copyright 2011 Mario Mulansky
 * Copyright 2012 Karsten Ahnert
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#include <boost/timer.hpp>
#include <boost/array.hpp>

#include <boost/numeric/odeint/stepper/runge_kutta4_classic.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/algebra/array_algebra.hpp>

#include "lorenz.hpp"

typedef boost::timer timer_type;

typedef boost::array< double , 3 > state_type;

using namespace boost::numeric::odeint;
//typedef boost::numeric::odeint::runge_kutta4_classic< state_type > rk4_odeint_type;
typedef runge_kutta4_classic< state_type , double , state_type , double ,
                              range_algebra, default_operations, never_resizer > rk4_odeint_type;


const int loops = 21;
const int num_of_steps = 20000000;
const double dt = 1E-10;


int main()
{
    double min_time = 1E6; // something big
    rk4_odeint_type stepper;
    for( int n=0; n<loops; n++ )
    {
        state_type x = {{ 8.5, 3.1, 1.2 }};
        double t = 0.0;
        timer_type timer;
        for( size_t i = 0 ; i < num_of_steps ; ++i )
        {
            stepper.do_step( lorenz(), x, t, dt );
            t += dt;
        }
        min_time = std::min( timer.elapsed() , min_time );
        std::clog << timer.elapsed() << '\t' << x[0] << std::endl;
    }
    std::cout << "Minimal Runtime: " << min_time << std::endl;
}
