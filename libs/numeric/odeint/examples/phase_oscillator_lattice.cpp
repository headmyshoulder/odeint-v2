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
#include <utility>

#include <boost/numeric/odeint.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace boost::numeric::odeint;

//[ phase_oscillator_lattice_system_function
typedef vector< double > container_type;


pair< double , double > calc_mean_field( const container_type &x )
{
	size_t n = x.size();
	double cos_sum = 0.0 , sin_sum = 0.0;
	for( size_t i=0 ; i<n ; ++i )
	{
		cos_sum += cos( x[i] );
		sin_sum += sin( x[i] );
	}
	cos_sum /= double( n );
	sin_sum /= double( n );

	double K = sqrt( cos_sum * cos_sum + sin_sum * sin_sum );
	double Theta = atan2( sin_sum , cos_sum );

	return make_pair( K , Theta );
}


struct phase_lattice
{
    container_type m_omega;
    double m_coupling;

    phase_lattice( const size_t n , double g = 1.0 ) : m_omega( n , 0.0 )
    {
        create_frequencies( g );
    }

    void create_frequencies( double g )
    {
    	boost::mt19937 rng;
    	boost::cauchy_distribution< double > cauchy( 0.0 , 0.001 );
    	for( size_t i=0 ; i<m_omega.size() ; ++i )
    	{
    		m_omega[i] = cauchy( rng );
    		cout << m_omega[i] << endl;
    	}
    	exit( -1 );
    }

    void set_coupling( double coupling ) { m_coupling = coupling; }

    double get_coupling( void ) const { return m_coupling; }

    void operator()( const container_type &x , container_type &dxdt , double t ) const
    {
    	pair< double , double > mean = calc_mean_field( x );

    	for( size_t i=0 ; i<x.size() ; ++i )
    		dxdt[i] = m_omega[i] + m_coupling * mean.first * sin( mean.second - x[i] );
    }
};
//]



//[ phase_oscillator_lattice_observer
struct statistics_observer
{
	double m_K_mean;
	size_t m_count;

	statistics_observer( void )
	: m_K_mean( 0.0 ) , m_count( 0 ) { }

    template< class State >
    void operator()( const State &x , double t )
    {
    	pair< double , double > mean = calc_mean_field( x );
    	m_K_mean += mean.first;
    	++m_count;
//    	cout << mean.first << "\t" << mean.second << "\t" << m_count << "\t" << m_K_mean << "\n";
    }

    double get_K_mean( void ) const { return ( m_count != 0 ) ? m_K_mean / double( m_count ) : 0.0 ; }

    void reset( void ) { m_K_mean = 0.0; m_count = 0; }
};
//]








int main( int argc , char **argv )
{
    //[ phase_oscillator_lattice_integration
//    const size_t n = 16384;
	const size_t n = 16;
    container_type x( n );
    generate( x.begin() , x.end() , drand48 );


    const double dt = 0.1;

    // gamma = 1
    // phase transition if coupling = 2
    phase_lattice lattice( n , 1.0 );
    statistics_observer obs;

    for( double coupling = 0.0 ; coupling < 5.0 ; coupling += 0.1 )
    {
    	lattice.set_coupling( coupling );
    	obs.reset();
    	integrate_const( explicit_rk4< container_type >() , boost::ref( lattice ) , x , 0.0 , 10.0 , dt , boost::ref( obs ) );
    	cout << coupling << "\t" << obs.get_K_mean() << "\n";
    }


    //]

    return 0;
}
