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

#include <boost/numeric/odeint/integrate/controller/conditional_stop.hpp>
#include <boost/numeric/odeint/integrate/controller/approximate.hpp>
#include <boost/numeric/odeint/integrate/integrate_conditional.hpp>

#include <boost/numeric/odeint.hpp>

#include <iostream>
#include <array>

namespace odeint = boost::numeric::odeint;

typedef std::array< double , 3 > state_type;


struct ball
{
    template< class State , class Time > 
    void operator()( const State &x , State &dxdt , Time t ) const
    {
        const double acceleration = 9.81;
        const double friction = 0.01;
        dxdt[0] = x[1];
        dxdt[1] = -acceleration - friction * x[1] ;
    }
};

struct ball_controller
{
    template< class Stepper , class Sys , class State , class Time >
    void init( Stepper &stepper , Sys sys , const State &s , Time t , Time dt )
    {
    }

    template< class State , class Time >
    bool stop( const State &s , Time &t )
    {
        return ( t > 100.0 );
    }

    template< class Stepper , class Sys , class State , class Time >
    void do_step( Stepper &stepper , Sys sys , State &x , Time &t , Time &dt )
    {
        stepper.do_step( sys , x , t , dt );
        t += dt;
    }

    template< class Stepper , class Sys , class State , class Time >
    void exit( Stepper &stepper , Sys sys , State &x , Time t , Time dt )
    {
    }
};

int main( int argc , char *argv[] )
{
    using namespace odeint;
    using namespace std;

    cout.precision( 14 );

    auto observer = []( const state_type &x , double t ) { cout << t << " " << x[0] << " " << x[1] << " " << "\n"; };

    // stepper concept
    if( false )
    {
        runge_kutta4< state_type > rk4;
        state_type x = {{ 10.0 , 0.0 }};
        integrate_conditional( rk4 , ball() , x , 0.0 , 0.01 ,
                               ball_controller() ,
                               observer );
    }




    return 0;
}
