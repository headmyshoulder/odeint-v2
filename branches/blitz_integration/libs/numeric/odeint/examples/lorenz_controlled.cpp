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
#include <fstream>
#include <vector>
#include <list>
#include <tr1/array>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/controlled_stepper_bs.hpp>

#define tab "\t"

using namespace std;
using namespace std::tr1;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

typedef array<double, 3> state_type;

size_t function_calls = 0;

void lorenz( state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
    function_calls++;
}


class output_observer
{

    ofstream m_file;

public:
    output_observer( const char* file_name )
        : m_file( file_name )
    { 
        m_file.precision(10);
        //m_file.setf(ios::fixed,ios::floatfield);
    }

    ~output_observer()
    {
        m_file.close();
    }

    void operator()( double t, state_type &x, ... )
    {
        m_file << t << tab << x[0] << tab << x[1] << tab << x[2] << endl;
    }
};

int main( int argc , char **argv )
{
    const double end_time = 25.0;

    const double eps_abs = 1E-10;
    const double eps_rel = 1E-10;;

    state_type x;
    x[0] = 1.0;
    x[1] = 0.0;
    x[2] = 20.0;

    stepper_rk5_ck< state_type > rk5;
    controlled_stepper_standard< stepper_rk5_ck< state_type > > controlled_rk5( rk5, eps_abs , eps_rel, 1.0, 1.0 );
    output_observer rk5_obs("lorenz_rk5.dat");
    size_t steps = integrate_adaptive( controlled_rk5, lorenz, x, 0.0, end_time, 1E-2, rk5_obs );

    clog << "RK5: " << steps << " steps. (" << function_calls << " function calls)" << endl;

    x[0] = 1.0;
    x[1] = 0.0;
    x[2] = 20.0;

    function_calls = 0;

    controlled_stepper_bs< state_type > controlled_bs(eps_abs, eps_rel, 1.0, 1.0);
    
    output_observer bs_obs("lorenz_bs.dat");
    steps = integrate_adaptive( controlled_bs, lorenz, x, 0.0, end_time, 1E-2, bs_obs );

    clog << "BS: " << steps << " steps. (" << function_calls << " function calls)" << endl;


    
    return 0;
}

/*
  Compile with
  g++ -Wall -O3 -I$BOOST_ROOT -I../../../../ lorenz_controlled.cpp
*/
