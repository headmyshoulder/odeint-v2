/* Boost numeric/odeint/examples/lorenz_array.cpp
 
 Copyright 2009 Karsten Ahnert

 Shows, the usage of odeint, and integrates the Lorenz equations,
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

const double dt = 0.01;
const size_t olen = 10000;

typedef array<double, 3> state_type;

void lorenz( state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

int main( int argc , char **argv )
{
    state_type x , x2;
    x[0] = 1.0;
    x[1] = 0.0;
    x[2] = 20.0;
    x2 = x;

    ode_step_euler< state_type > euler;
    ode_step_runge_kutta_4< state_type > rk4;

    double t = 0.0;
    for( size_t oi=0 ; oi<olen ; ++oi,t+=dt )
    {
	cout << t << tab << x[0] << tab << x[1] << tab << x[2] << tab;
	cout << x2[0] << tab << x2[1] << tab << x2[2] << endl;
	euler.next_step( lorenz , x , t , dt );
	rk4.next_step( lorenz , x2 , t , dt );
    }

    return 0;
}



/*
  Compile with
  g++ -Wall -I$BOOST_ROOT -I../../../../ lorenz_array.cpp
*/
