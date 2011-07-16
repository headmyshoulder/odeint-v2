/*
 * two_dimensional_phase_lattice.cpp
 *
 * This example show how one can use matrices as state types in odeint.
 *
 *  Created on: Jul 15, 2011
 *      Author: karsten
 *
 * Copyright 2009 Karsten Ahnert and Mario Mulansky.
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 */

#include <iostream>

#include <boost/numeric/odeint.hpp>


using namespace boost::numeric::odeint;

typedef boost::numeric::ublas::matrix< double , boost::numeric::ublas::basic_row_major<> , boost::numeric::ublas::unbounded_array< double > > state_type;

struct two_dimensional_phase_lattice
{
    void operator()( const state_type &x , state_type &dxdt , double t ) const
    {
    }
};

int main( int argc , char **argv )
{
    state_type x( 512 , 512 , 0.0 );

    integrate_const( explicit_rk4< state_type >() , two_dimensional_phase_lattice() ,
            x , 0.0 , 10000.0 , 0.1 );


    return 0;
}
