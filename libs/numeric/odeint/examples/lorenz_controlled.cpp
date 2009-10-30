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

const size_t olen = 10000;

const double eps_abs = 1E-3;
const double eps_rel = 1E-3;

const double min_dt = 1E-10;

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

    ode_step_euler< state_type > euler;
    step_controller_standard< state_type, double > controlled_euler(eps_abs, 
								    eps_rel, 
								    1.0, 1.0);
    
    double t = 0.0;
    double dt = 0.01;
    controlled_step_result result;

    cout.precision(5);
    cout.setf(ios::fixed,ios::floatfield);
    

    while( (t<10.0) && (dt > min_dt) ) {
	//cout << "current state: " ;
	cout << t << tab << dt << tab << x[0] << tab << x[1] << tab << x[2] << endl;
	result = controlled_euler.controlled_step( euler , lorenz , x , t , dt );
	while( result != SUCCESS ) {
	    result = controlled_euler.controlled_step( euler , lorenz , x , t , dt );
	    if( dt < min_dt ) break;
	}
	//cout<<"SUCCESS with dt="<<dt<<endl;
    }

    return 0;
}

/*
  Compile with
  g++ -Wall -I$BOOST_ROOT -I../../../../ lorenz_array.cpp
*/
