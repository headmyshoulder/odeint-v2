#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <boost/numeric/odeint.hpp>

#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;


typedef std::tr1::array< double , 2 > state_type;


const double gam = .15;

void harmonic_oscillator(const state_type &x, state_type &dxdt, const double t)
{
    dxdt[0] = x[1];
    dxdt[1] = -x[0] - gam*x[1];
}


int main( int argc , char **argv )
{
    const double dt = 0.01;
    const size_t olen = 10000;
    state_type x1 = {{ 1.0, 0.0 }};
    state_type x2 = {{ 1.0, 0.0 }};
    state_type x3 = {{ 1.0, 0.0 }};
    state_type x4 = {{ 1.0, 0.0 }};
    state_type x5 = {{ 1.0, 0.0 }};

    stepper_euler< state_type > stepper_euler;
    stepper_midpoint< state_type > stepper_mp;
    stepper_rk4< state_type > stepper_rk4;

    double t = 0.0;
    for( size_t oi=0 ; oi<olen ; ++oi,t+=dt )
    {
        stepper_euler.next_step( harmonic_oscillator , x1 , t , dt );
        stepper_mp.next_step( harmonic_oscillator , x2 , t , dt , 2 );
        stepper_mp.next_step( harmonic_oscillator , x3 , t , dt , 4 );
        stepper_mp.next_step( harmonic_oscillator , x4 , t , dt , 8 );
        stepper_rk4.next_step( harmonic_oscillator , x5 , t , dt );
        cout<< t << tab << x1[0]*x1[0] + x1[1]*x1[1];
        cout<< tab << x2[0]*x2[0] + x2[1]*x2[1];
        cout<< tab << x3[0]*x3[0] + x3[1]*x3[1];
        cout<< tab << x4[0]*x4[0] + x4[1]*x4[1];
        cout<< tab << x5[0]*x5[0] + x5[1]*x5[1] << endl;
    }
}

/*
  Compile with
  g++ -Wall -O3 -I$BOOST_ROOT -I../../../../../ stepper_midpoint.cpp
*/
