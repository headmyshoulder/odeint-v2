/* Boost numeric/odeint/examples/lorenz_controlled.cpp
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Shows the usage of odeint by integrating the Lorenz equations,
 a deterministic chaotic system

 dx/dt = sigma * ( x - y)
 dy/dt = R*x - y - x*z
 dz/dt = x*y - b*z

 mit sigma = 10, r=28, b = 8/3

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <iostream>
#include <vector>
#include <list>
#include <tr1/array>

#include <boost/numeric/odeint.hpp>

#define tab "\t"

using namespace std;
using namespace std::tr1;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

const size_t olen = 10000;

const double eps_abs = 1E-6;
const double eps_rel = 1E-7;

typedef array<double, 3> state_type;

void lorenz( state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

void print_state( double t, state_type &x, ... )
{
    cout << t << tab << x[0] << tab << x[1] << tab << x[2] << endl;
}

int main( int argc , char **argv )
{
    state_type x;
    x[0] = 1.0;
    x[1] = 0.0;
    x[2] = 20.0;

    ode_step_half_step< ode_step_euler< state_type > > euler;
    step_controller_standard< state_type, double > controller( eps_abs , eps_rel, 1.0, 1.0);
    
    cout.precision(5);
    cout.setf(ios::fixed,ios::floatfield);
    
    size_t steps = integrate( euler, lorenz, controller, 0.0, 1E-4, x, 10.0, print_state );

    clog << "Number of steps: " << steps << endl;

    return 0;
}

/*
  Compile with
  g++ -Wall -O3 -I$BOOST_ROOT -I../../../../ lorenz_controlled.cpp
*/
