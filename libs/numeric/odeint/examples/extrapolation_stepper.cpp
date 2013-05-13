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
    
typedef extrapolation_stepper< 4 , state_type > stepper_type4;
typedef extrapolation_stepper< 6 , state_type > stepper_type6;
typedef extrapolation_stepper< 8 , state_type > stepper_type8;

int main()
{
    double dt = 0.5;
    double t = 0.0;
    while( dt > 1E-4 )
    {
        state_type x1 = {{ 0.0 , 1.0 }};
        state_type err1;
        stepper_type4 stepper4;
        stepper4.do_step( osc() , x1 , t , dt , err1 );

        state_type x2 = {{ 0.0 , 1.0 }};
        state_type err2;
        stepper_type6 stepper6;
        stepper6.do_step( osc() , x2 , t , dt , err2 );

        state_type x3 = {{ 0.0 , 1.0 }};
        state_type err3;
        stepper_type8 stepper8;
        stepper8.do_step( osc() , x3 , t , dt , err3 );

        using std::abs;
        std::cout << dt << '\t' << abs(x1[0] - sin(dt)) << '\t' << abs( err1[0] ) << '\t';
        std::cout << abs(x2[0] - sin(dt)) << '\t' << abs( err2[0] ) << '\t';
        std::cout << abs(x3[0] - sin(dt)) << '\t' << abs( err3[0] ) << std::endl;
        
        dt /= 1.257;
    }
}
