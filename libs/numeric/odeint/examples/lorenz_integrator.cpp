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
#include <iterator>
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

const double eps_abs = 1E-6;
const double eps_rel = 1E-7;

const size_t time_points = 101;

typedef array<double, 3> state_type;

void lorenz( state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

int main( int argc , char **argv )
{
    state_type x1;
    x1[0] = 1.0;
    x1[1] = 0.0;
    x1[2] = 20.0;
    state_type x2(x1);
    state_type x3(x1);


    vector<state_type> x1_t_vec;
    vector<double> t1_vec;
    vector<state_type> x2_t_vec;
    vector<double> t2_vec;
    vector<state_type> x3_t_vec;
    vector<double> t3_vec;

    stepper_half_step< stepper_euler< state_type > > euler;
    controlled_stepper_standard< stepper_half_step< stepper_euler< state_type > > >
        euler_controlled( euler , eps_abs, eps_rel, 1.0, 1.0);
    size_t steps = integrate( euler_controlled, lorenz, x1, 
                              0.0, 10.0, 1E-4, 
                              back_inserter(t1_vec),
                              back_inserter(x1_t_vec));

    clog << "Euler Half Step: #steps " << steps << endl;

    stepper_half_step< stepper_rk4< state_type > > rk4;
    controlled_stepper_standard< stepper_half_step< stepper_rk4< state_type > > >
        rk4_controlled( rk4 , eps_abs, eps_rel, 1.0, 1.0);
    steps = integrate( rk4_controlled, lorenz, x2, 0.0, 10.0, 1E-4, 
                       back_inserter(t2_vec),
                       back_inserter(x2_t_vec));

    clog << "RK4 Half Step: #steps " << steps << endl;


    stepper_rk5_ck< state_type > rk5;
    controlled_stepper_standard< stepper_rk5_ck< state_type > > 
        rk5_controlled( rk5 , eps_abs, eps_rel, 1.0, 1.0);
    steps = integrate( rk5_controlled, lorenz, x3, 0.0, 10.0, 1E-4, 
                       back_inserter(t3_vec),
                       back_inserter(x3_t_vec));
    
    clog << "RK5 Cash-Karp: #steps: " << steps << endl;

    cout.precision(16);
    cout.setf(ios::fixed,ios::floatfield);

    cout << t1_vec.size() << tab << t2_vec.size() << tab << t3_vec.size() << endl;
    

    //cout << "current state: " ;
    cout << (x1_t_vec.back())[0] << tab << (x1_t_vec.back())[1] << tab << (x1_t_vec.back())[2] << tab;
    cout << x2_t_vec.back()[0] << tab << x2_t_vec.back()[1] << tab << x2_t_vec.back()[2] << tab;
    cout << x3_t_vec.back()[0] << tab << x3_t_vec.back()[1] << tab << x3_t_vec.back()[2] << endl;

    return 0;
}

/*
  Compile with
  g++ -Wall -I$BOOST_ROOT -I../../../../ lorenz_integrator.cpp
*/
