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
    state_type x1 = {{ 1.0, 0.0 }};

    const size_t step_count = 100;
    controlled_stepper_bs< state_type > stepper_bs(0.0, 1E-11, 1.0, 1.0);
    double t = 0.0;
    double dt = 0.01;
    cout << " Everything initialized - starting time evolution " << endl;
    cout.precision(16);
    clog.precision(16);
    for( size_t step=0; step < step_count; step++ )
    {
        stepper_bs.try_step(harmonic_oscillator, x1, t, dt);
        //clog << " ####################################################### " << endl;
        cout << t << tab << dt << tab << x1[0] << tab << x1[1] << endl;
        //clog << " ####################################################### " << endl;
    }
}

/*
  Compile with
  g++ -Wall -O3 -I$BOOST_ROOT -I../../../../../ stepper_bs.cpp
*/

