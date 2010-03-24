/* Boost numeric/odeint/examples/coupled_vdp.cpp
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Shows the usage of odeint by integrating the equations of a 
 pendulum with horizontally vibrating pivot:

 d^2x/dt^2 + sin(x) + alpha*x = a*omega^2*sin(omega*t)*cos(x)

 for large enough omega >sim 20 two new fixpoints (of the 
 slow dynamics) arise, that can be seen in the simulations as well

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

#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;


typedef std::tr1::array< double , 2 > state_type;

const double alpha = 0.1;
const double omega = 20;
const double a = 0.1;

/* 
   Defines the right hand side f(x,t) of the dynamical equations dx/dt = f(x,t) 
   x consists of x=(x, dx/dt) and f has explicit time dependence
*/
void my_system( state_type &x , state_type &dxdt , double t ) 
{
    dxdt[0] = x[1];
    dxdt[1] = -sin(x[0]) - alpha*x[1] + a*omega*omega*sin(omega*t)*cos(x[0]);
}

int main( int argc , char **argv )
{
    state_type x = {{ 1.0, 0.0 }};

    vector<double> times;
    vector<state_type> x_t_vec;
    
    stepper_half_step< stepper_rk4< state_type > > stepper;

    controlled_stepper_standard
        < stepper_half_step< stepper_rk4< state_type > >
        > controlled_stepper( 1E-6 , 1E-7 , 1.0 , 1.0 );

    size_t steps = integrate( controlled_stepper, my_system, x, 
                              0.0, 100.0,1E-4, 
                              back_inserter(times),
                              back_inserter(x_t_vec));

    clog << "Steps: " << steps << endl;

    cout.precision(5);
    cout.setf(ios::fixed,ios::floatfield);
    

    for( size_t i=0; i<times.size(); i++ )
    {
        //cout << "current state: " ;
        cout << times[i] << tab;
        cout << x_t_vec[i][0] << tab << x_t_vec[i][1] << endl;
    }

    return 0;

}
