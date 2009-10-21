/* Boost numeric/odeint/examples/lorenz.cpp
 
 Copyright 2009 Karsten Ahnert

 Shows, the usage of odeint, and integrates the Lorenz equations,
 a deterministic chaotic system. In this version, the system state
 is represented by a c-style double[3] array. The Lorenz system is
 defined as

 dx/dt = sigma * ( x - y)
 dy/dt = R*x - y - x*z
 dz/dt = x*y - b*z

 with sigma = 10, r=28, b = 8/3

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <iostream>

#include <boost/numeric/odeint.hpp>

#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

const double dt = 0.01;
const size_t olen = 10000;

void lorenz( double* x , double* dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

int main( int argc , char **argv )
{
    double x[3];
    x[0] = 1.0;
    x[1] = 0.0;
    x[2] = 0.0;

    ode_step_euler< double* > euler;

    double t = 0.0;
    for( size_t oi=0 ; oi<olen ; ++oi,t+=dt )
    {
	cout << t << tab << x[0] << tab << x[1] << tab << x[2] << endl;
	euler.next_step( lorenz , x , 3, t , dt );
    }

    return 0;
}

/*
  Compile with
  g++ -Wall -I../../../../ lorenz_carray.cpp
*/
