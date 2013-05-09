/* Boost libs/numeric/odeint/examples/extrapolation_stepper.cpp

 Copyright 2009-2013 Karsten Ahnert
 Copyright 2009-2013 Mario Mulansky

 example usage of the extrapolation stepper with arbirary order

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#include <iostream>
#include <cmath>

#include <boost/numeric/odeint.hpp>
#include <boost/ref.hpp>

using namespace boost::numeric::odeint;


typedef boost::array< double , 2 > state_type;

// harmonic oscillator, analytic solution x[0] = sin( t )
struct osc
{
    void operator()( const state_type &x , state_type &dxdt , const double t ) const
    {
        dxdt[0] = x[1];
        dxdt[1] = -x[0];
    }
};
    
typedef extrapolation_stepper< 3 , state_type > stepper_type3;
typedef extrapolation_stepper< 5 , state_type > stepper_type5;
typedef extrapolation_stepper< 7 , state_type > stepper_type7;

int main()
{
    double dt = 0.5;
    double t = 0.0;
    while( dt > 1E-4 )
    {
        state_type x_out;
        state_type x3 = {{ 0.0 , 1.0 }};
        stepper_type3 stepper3;
        stepper3.do_step( osc() , x3 , t , x_out , dt );

        state_type x5 = {{ 0.0 , 1.0 }};
        stepper_type5 stepper5;
        stepper5.do_step( osc() , x5 , t , dt );

        state_type x7 = {{ 0.0 , 1.0 }};
        stepper_type7 stepper7;
        stepper7.do_step( osc() , x7 , t , dt );

        std::cout << dt << '\t' << std::abs(x_out[0] - sin(dt)) << '\t';
        std::cout << std::abs(x5[0] - sin(dt)) << '\t';
        std::cout << std::abs(x7[0] - sin(dt)) << std::endl;
        
        dt /= 1.257;
    }
}
