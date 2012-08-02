/*
 * odeint_rk4_lorenz.cpp
 *
 *  Created on: May 11, 2011
 *      Author: mario
 */

#include <algorithm>
#include <array>

#include <boost/range/algorithm.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/numeric.hpp>
#include <boost/timer.hpp>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/generation.hpp>
#include <boost/numeric/odeint/iterator/const_step_iterator.hpp>
#include <boost/numeric/odeint/iterator/const_step_time_iterator.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef std::array< double , 3 > state_type;

std::ostream& operator<<( std::ostream &out , const state_type &x )
{
    out << x[0] << "\t" << x[1] << "\t" << x[2];
    return out;
}

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

struct lorenz
{
    template< class State , class Deriv >
    void operator()( const State &x , Deriv &dxdt , double t ) const
    {
        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = -b * x[2] + x[0] * x[1];
    }
};

state_type test1_iterator( double t_end = 1000000.0 , double dt = 0.01 )
{
    runge_kutta4< state_type > stepper;
    state_type x = {{ 10.0 , 10.0 , 10.0 }};
    auto first = make_const_step_time_iterator_begin( stepper , lorenz() , x , 0.0 , dt );
    auto last = make_const_step_time_iterator_end( stepper , lorenz() , x , t_end , dt );
    for( ; first != last ; )
        ++first;
    return x;
}

state_type test1_integrate( double t_end = 1000000.0 , double dt = 0.01 )
{
    runge_kutta4< state_type > stepper;
    state_type x = {{ 10.0 , 10.0 , 10.0 }};
    integrate_const( stepper , lorenz() , x , 0.0 , t_end , dt );
    return x;
}



int main( int argc , char **argv )
{
    boost::timer timer;

    cout << "Test 1 : " << endl;
    timer.restart();
    state_type x = test1_iterator( 1000000.0 );
    cout << "\tIterator : " << x << " elapsed time : " << timer.elapsed() << " s" << endl;

    timer.restart();
    x = test1_integrate( 1000000.0 );
    cout << "\tIntegrate : " << x << " elapsed time : " << timer.elapsed() << " s" << endl;
    
    return 0;
}
