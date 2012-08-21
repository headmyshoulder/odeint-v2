/*
 * odeint_rk4_lorenz.cpp
 *
 * Copyright 2009-2012 Karsten Ahnert
 * Copyright 2009-2012 Mario Mulansky
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
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
typedef std::vector< double  > state_type2;

std::ostream& operator<<( std::ostream &out , const state_type &x )
{
    out << x[0] << "\t" << x[1] << "\t" << x[2];
    return out;
}

std::ostream& operator<<( std::ostream &out , const state_type2 &x )
{
    if( !x.empty() ) out << x[0];
    for( size_t i=1 ; i<x.size() ; ++i ) out << "\t" << x[i];
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

boost::timer timer;

void test1_iterator( double t_end = 1000000.0 , double dt = 0.01 )
{
    cout << "Test 1 Iterator with array : ";
    timer.restart();

    runge_kutta4< state_type > stepper;
    state_type x = {{ 10.0 , 10.0 , 10.0 }};
    auto first = make_const_step_iterator_begin( stepper , lorenz() , x , 0.0 , t_end , dt );
    auto last = make_const_step_iterator_end( stepper , lorenz() , x );
    for( ; first != last ; )
        ++first;

    cout << "\t" << x << " elapsed time : " << timer.elapsed() << " s" << endl;
}

void test1_time_iterator( double t_end = 1000000.0 , double dt = 0.01 )
{
    cout << "Test 1 Time iterator with array : ";
    timer.restart();

    runge_kutta4< state_type > stepper;
    state_type x = {{ 10.0 , 10.0 , 10.0 }};
    auto first = make_const_step_time_iterator_begin( stepper , lorenz() , x , 0.0 , t_end , dt );
    auto last = make_const_step_time_iterator_end( stepper , lorenz() , x );
    for( ; first != last ; )
        ++first;

    cout << "\t" << x << " elapsed time : " << timer.elapsed() << " s" << endl;
}

void test1_integrate( double t_end = 1000000.0 , double dt = 0.01 )
{
    cout << "Test 1 Integrate with array : ";
    timer.restart();

    runge_kutta4< state_type > stepper;
    state_type x = {{ 10.0 , 10.0 , 10.0 }};
    integrate_const( stepper , lorenz() , x , 0.0 , t_end , dt );

    cout << "\t" << x << " elapsed time : " << timer.elapsed() << " s" << endl;
}




void test2_iterator( double t_end = 1000000.0 , double dt = 0.01 )
{
    cout << "Test 2 Iterator with vector : ";
    timer.restart();

    runge_kutta4< state_type2 > stepper;
    state_type2 x = { 10.0 , 10.0 , 10.0 };
    auto first = make_const_step_iterator_begin( stepper , lorenz() , x , 0.0 , t_end , dt );
    auto last = make_const_step_iterator_end( stepper , lorenz() , x );
    for( ; first != last ; )
        ++first;

    cout << "\t" << x << " elapsed time : " << timer.elapsed() << " s" << endl;
}

void test2_time_iterator( double t_end = 1000000.0 , double dt = 0.01 )
{
    cout << "Test 2 Time iterator with vector : ";
    timer.restart();

    runge_kutta4< state_type2 > stepper;
    state_type2 x = { 10.0 , 10.0 , 10.0 };
    auto first = make_const_step_time_iterator_begin( stepper , lorenz() , x , 0.0 , t_end , dt );
    auto last = make_const_step_time_iterator_end( stepper , lorenz() , x );
    for( ; first != last ; )
        ++first;

    cout << "\t" << x << " elapsed time : " << timer.elapsed() << " s" << endl;
}

void test2_integrate( double t_end = 1000000.0 , double dt = 0.01 )
{
    cout << "Test 2 Integrate with vector : ";
    timer.restart();
    runge_kutta4< state_type2 > stepper;
    state_type2 x = { 10.0 , 10.0 , 10.0 };
    integrate_const( stepper , lorenz() , x , 0.0 , t_end , dt );
    cout << "\t" << x << " elapsed time : " << timer.elapsed() << " s" << endl;
}



int main( int argc , char **argv )
{
    test1_iterator( 1000000.0 );
    test1_time_iterator( 1000000.0 );
    test1_integrate( 1000000.0 );

    test2_iterator( 1000000.0 );
    test2_time_iterator( 1000000.0 );
    test2_integrate( 1000000.0 );


    
    return 0;
}
