/* Boost numeric/odeint/examples/lorenz_integrator.cpp
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Shows the usage of odeint integrator by integrating the Lorenz equations,
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

const double eps_abs = 1E-3;
const double eps_rel = 1E-3;

const size_t time_points = 100;

typedef array<double, 3> state_type;

void lorenz( state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

int main( int argc , char **argv )
{
    state_type x;
    x[0] = 1.0;
    x[1] = 0.0;
    x[2] = 20.0;

    vector<state_type> x_t_vec(time_points);
    vector<double> times(time_points);
    for( size_t i=0; i<time_points; i++ ) {
	times[i] = 0.1*i;
    }

    ode_step_euler< state_type > euler;
    size_t steps = integrate( euler, lorenz, x, times, x_t_vec);

    clog << "Steps: " << steps << endl;

    cout.precision(5);
    cout.setf(ios::fixed,ios::floatfield);
    

    for( size_t i=0; i<time_points; i++ ) {
	//cout << "current state: " ;
	cout << times[i] << tab;
	cout << x_t_vec[i][0] << tab << x_t_vec[i][1] << tab << x_t_vec[i][2] << endl;
    }

    return 0;
}

/*
  Compile with
  g++ -Wall -I$BOOST_ROOT -I../../../../ lorenz_integrator.cpp
*/
