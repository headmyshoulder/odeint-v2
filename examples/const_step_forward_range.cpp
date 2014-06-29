/*
 * const_step_iterator.cpp
 *
 * Copyright 2012-2013 Karsten Ahnert
 * Copyright 2013 Mario Mulansky
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * several examples for using iterators
 */


#include <iostream>
#include <iterator>
#include <utility>
#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>

#include <boost/range/algorithm.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/numeric.hpp>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/generation.hpp>
#include <boost/numeric/odeint/range/const_step_forward_range.hpp>


#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;

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



int main( int argc , char **argv )
{
    typedef std::array< double , 3 > state_type;

    runge_kutta4< state_type > stepper;
    state_type x = {{ 10.0 , 10.0 , 10.0 }};

    auto r = range::make_const_step_forward_range( stepper , x , lorenz() , 0.0 , 0.1 , 0.01 );
    auto first = r.begin();
    auto last = r.end();

    r.print_state();
    std::cout << "\n\n\n\n";

    // size_t count = 0;
    // std::ofstream fout { "lorenz.dat" };
    // for( ; first != last ; ++count )
    // {
    //     state_type const& x = *first++;
    //     std::cout << "\n\n\n\n";
    //     fout << x[0] << tab << x[1] << tab << x[2] << "\n";

    //     // if( count == 3 ) break;
    // }

    boost::range::for_each( r ,
                            []( const state_type &x ) {
                                std::cout << x[0] << tab << x[1] << tab << x[2] << "\n"; } );

    return 0;
}
