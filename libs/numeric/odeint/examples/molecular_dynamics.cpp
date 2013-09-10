/*
 [auto_generated]
 libs/numeric/odeint/examples/molecular_dynamics.cpp

 [begin_description]
 tba.
 [end_description]

 Copyright 2009-2012 Karsten Ahnert
 Copyright 2009-2012 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/numeric/odeint.hpp>

#include <vector>
#include <iostream>
#include <random>

using namespace boost::numeric::odeint;

struct md_system
{
    static const size_t n = 128;
    typedef std::vector< double > state_type;
    
    md_system( double a = 1.0 , double gamma = 0.01 ) : m_a( a ) , m_gamma( gamma ) { }
    
    static void init_state_type( state_type &x ) { x.resize( 2 * n ); }
    
    void operator()( state_type const& x , state_type const& v , state_type &a , double t ) const
    {
        for( size_t i=0 ; i<n ; ++i )
        {
            double r2 = x[i] * x[i] + x[n+1] * x[n+i];
            double r = std::sqrt( r2 );
            a[     i ] = - m_a * r - m_gamma * v[     i ] ;
            a[ n + i ] = - m_a * r - m_gamma * v[ n + i ] ;
        }
    }
    
    double m_a;
    double m_gamma;
};

int main( int argc , char *argv[] )
{
    const size_t n = md_system::n;
    typedef md_system::state_type state_type;
    
    
    std::mt19937 rng;
    std::normal_distribution<> x_dist( 0.0 , 1.0 );
    
    state_type x , v;
    md_system::init_state_type( x );
    md_system::init_state_type( v );
    
    
    for( size_t i=0 ; i<n ; ++i )
    {
        x[  i] = x_dist( rng );
        x[n+i] = x_dist( rng );
        v[  i] = 0.0;
        v[n+i] = 0.0;
    }
    
    for( size_t i=0 ; i<10000 ; ++i )
    {
    }
    
    integrate_const( velocity_verlet< state_type >() , md_system() , std::make_pair( std::ref( x ) , std::ref( v ) ) , 0.0 , 1.0 , 0.01 ,
        []( 
    );
    
    return 0;
}
