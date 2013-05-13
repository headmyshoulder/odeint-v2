/*
 * Copyright 2009-2013 Karsten Ahnert
 * Copyright 2009-2013 Mario Mulansky
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * Example comparing performance of th extrapolation stepper and the rkdopri5
*/

#include <iostream>
#include <cmath>

#include <boost/array.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/timer.hpp>

#include <boost/numeric/odeint.hpp>

typedef boost::array< double , 3 > state_type;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

void lorenz( const state_type &x , state_type &dxdt , const double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = -b * x[2] + x[0] * x[1];
}

using namespace boost::numeric::odeint;
using namespace boost::accumulators;
using namespace std;

typedef accumulator_set<
    double , stats< tag::mean , tag::variance >
    > accumulator_type;

typedef runge_kutta_dopri5<state_type> dopri5_stepper;
typedef extrapolation_stepper< 6 , state_type > extrapol_stepper;

ostream& operator<<( ostream& out , accumulator_type &acc )
{
    out << boost::accumulators::mean( acc ) << '\t';
//    out << boost::accumulators::variance( acc ) << tab;
    return out;
}

typedef boost::timer timer_type;

int main()
{
    const size_t loops = 20;

    accumulator_type acc;
    timer_type timer;
    int steps;

    cout << "Dopri5: " << endl;

    for( size_t n=0 ; n<loops ; ++n )
    {
        timer.restart();
        state_type x = {{ 10.0 , 5.0 , 5.0 }};
        steps = integrate_adaptive( make_controlled<dopri5_stepper>( 1E-12 , 1E-12 ) , lorenz , x ,
                                    0.0 , 500.0 , 0.1 );
        if( n>0 )
        {   // take first run as transient
            acc(timer.elapsed());
            cout << timer.elapsed() << '\t' << steps << endl;
        }
    }
    double time_dopri5 = boost::accumulators::mean(acc);

    cout << "6th order Extrapolation: " << endl;
    acc = accumulator_type();
    for( size_t n=0 ; n<loops ; ++n )
    {
        timer.restart();
        state_type x = {{ 10.0 , 5.0 , 5.0 }};
        steps = integrate_adaptive( make_controlled<extrapol_stepper>( 1E-12 , 1E-12 ) , lorenz , x ,
                                    0.0 , 500.0 , 0.1 );
        if( n>0 )
        {   // take first run as transient
            acc(timer.elapsed());
            cout << timer.elapsed() << '\t' << steps << endl;
        }
    }
    double time_extra = boost::accumulators::mean(acc);
    cout <<" time for dopri5: " << time_dopri5 << " , time for extrapolation: " << time_extra << endl;
}
