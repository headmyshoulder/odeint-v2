/* Boost numeric/odeint/examples/lorenz_stepper.cpp
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Testing of the Hamiltonian solvers

 Furthermore, the usage of std::tr1::array and std::vector in odeint is
 shown and the performance of both containers is compared.

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
using namespace boost::numeric::odeint;


typedef std::tr1::array< double , 1 > state_type;

void harm_osc( state_type& q , state_type &dpdt )
{
    dpdt[0] = - q[0];
}




int main( int argc , char **argv )
{
    const double dt = 0.25;
    const size_t olen = 10000;
    
    state_type q1 = {{ 1.0 }} , p1 = {{ 0.0 }};
    state_type q2 = {{ 1.0 }} , p2 = {{ 0.0 }};

    hamiltonian_stepper_euler< state_type > euler;
    hamiltonian_stepper_rk< state_type > rk;

    double t = 0.0;
    for( size_t oi=0 ; oi<olen ; ++oi,t+=dt )
    {
        euler.do_step( harm_osc , q1 , p1 , dt );
        rk.do_step( harm_osc , q2 , p2 , dt );
	cout << t << tab << q1[0] << tab << p1[0] << tab;
	cout << q2[0] << tab << p2[0] << "\n";
    }

    return 0;
}



/*
  Compile with
  g++ -Wall -O3 -I$BOOST_ROOT -I../../../../ lorenz_stepper.cpp
*/
