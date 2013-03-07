/*
 [auto_generated]
 libs/numeric/odeint/examples/integrate_conditional.cpp

 [begin_description]
 tba.
 [end_description]

 Copyright 2009-2012 Karsten Ahnert
 Copyright 2009-2012 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/numeric/odeint/integrate/controller/adaptive_stop.hpp>
#include <boost/numeric/odeint/integrate/controller/adaptive_approximate.hpp>
#include <boost/numeric/odeint/integrate/integrate_conditional.hpp>

#include <boost/numeric/odeint.hpp>

#include <iostream>
#include <array>

namespace odeint = boost::numeric::odeint;

typedef std::array< double , 3 > state_type;


struct lorenz
{
    template< class State , class Time > 
    void operator()( const State &x , State &dxdt , Time t ) const
    {
        const double sigma = 10.0;
        const double R = 28.0;
        const double b = 8.0 / 3.0;
        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = -b * x[2] + x[0] * x[1];
    }
};

int main( int argc , char *argv[] )
{
    using namespace odeint;
    using namespace std;

    cout.precision( 14 );

    auto observer = []( const state_type &x , double t ) { cout << t << " " << x[0] << " " << x[1] << " " << x[2] << "\n"; };

    // stepper concept
    if( false )
    {
        runge_kutta4< state_type > rk4;
        state_type x = {{ 10.0 , 10.0 , 10.0 }};
        integrate_conditional( rk4 , lorenz() , x , 0.0 , 0.01 ,
                               make_adaptive_stop( []( const state_type &x , double t ) -> bool { return t > 0.995; } ) ,
                               observer );
    }

    // stepper concept
    {
        runge_kutta4< state_type > rk4;
        state_type x = {{ 10.0 , 10.0 , 10.0 }};
        integrate_conditional( rk4 , lorenz() , x , 0.0 , 0.01 ,
                               make_adaptive_stop( []( const state_type &x , double t ) -> bool { return x[0] < -10.0; } ) ,
                               observer );
    }


    // stepper concept
    {
        runge_kutta4< state_type > rk4;
        state_type x = {{ 10.0 , 10.0 , 10.0 }};
        integrate_conditional( rk4 , lorenz() , x , 0.0 , 0.01 ,
                               make_adaptive_approximate( []( const state_type &x , double t ) -> bool { return x[0] < -10.0; } , 1.0e-10 ) ,
                               observer );
    }




    return 0;
}
