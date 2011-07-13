/*
 * phase_oscillator_lattice.cpp
 *
 * Demonstrates the phase transition from an unsynchronized to an synchronized state.
 *
 *  Created on: Jul 13, 2011
 *      Author: karsten
 *
 * Copyright 2009 Karsten Ahnert and Mario Mulansky.
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 */

#include <iostream>
#include <tr1/array>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

//[ phase_oscillator_lattice_system_function
typedef vector< double > container_type;

struct phase_lattice
{
    container_type m_omega;
    double m_coupling;

    phase_lattice( const size_t n , double sigma = 1.0 ) : m_omega( n , 0.0 )
    {
        set_omega( sigma );
    }

    void set_omega( double sigma )
    {
        generate( m_omega.begin() , m_omega.end() , drand48 );
    }

    void set_coupling( double coupling ) { m_coupling = coupling; }

    double get_coupling( void ) const { return m_coupling; }

    void operator()( const container_type &x , container_type &dxdt , double t ) const
    {
    }
};
//]



//[ phase_oscillator_lattice_observer
struct statistics_observer
{
    template< class State >
    void operator()( const State &x , double t ) const
    {
    }
};
//]








int main( int argc , char **argv )
{
    //[ phase_oscillator_lattice_integration
    const size_t n = 16384;
    container_type x( n );
    generate( x.begin() , x.end() , drand48 );


    const double dt = 0.1;

    phase_lattice lattice( n );
    statistics_observer obs;

    integrate_const( explicit_rk4< container_type >() , boost::ref( lattice ) , x , 0.0 , 10000.0 , dt , obs );
    //]

    return 0;
}
