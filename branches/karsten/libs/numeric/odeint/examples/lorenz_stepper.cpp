/* Boost numeric/odeint/examples/lorenz_stepper.cpp
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Shows the usage of odeint, and integrates the Lorenz equations,
 a deterministic chaotic system

 dx/dt = sigma * ( x - y)
 dy/dt = R*x - y - x*z
 dz/dt = x*y - b*z

 with sigma = 10, r=28, b = 8/3

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


typedef std::vector< double > state_type1;
typedef std::tr1::array< double , 3 > state_type2;


const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

void lorenz1( state_type1 &x , state_type1 &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

void lorenz2( state_type2 &x , state_type2 &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}




int main( int argc , char **argv )
{
    const double dt = 0.01;
    const size_t olen = 100000000;
    
    state_type1 x1(3);
    x1[0] = 1.0;
    x1[1] = 0.0;
    x1[2] = 0.0;
    state_type2 x2 = {{ 1.0 , 0.0 , 0.0 }};

    stepper_rk4< state_type1 > stepper1( x1 );
    stepper_rk4< state_type2 > stepper2( x2 );

    clock_t start , end;
    double t;

    start= clock();
    t = 0.0;
    for( size_t oi=0 ; oi<olen ; ++oi,t+=dt )
        stepper1.do_step( lorenz1 , x1 , t , dt );
    end = clock();
    cout << "vector : " << double ( end - start ) / double( CLOCKS_PER_SEC ) << endl;
    cout << "x: "<<x1[0]<<tab<<x1[1]<<tab<<x1[2]<<endl;


    start= clock();
    t = 0.0;
    for( size_t oi=0 ; oi<olen ; ++oi,t+=dt )
        stepper2.do_step( lorenz2 , x2 , t , dt );
    end = clock();
    cout << "array : " << double ( end - start ) / double( CLOCKS_PER_SEC ) << endl;
    cout << "x: "<<x2[0]<<tab<<x2[1]<<tab<<x2[2]<<endl;

    return 0;
}



/*
  Compile with
  g++ -Wall -O3 -I$BOOST_ROOT -I../../../../ lorenz_stepper.cpp
*/
