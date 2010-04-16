/* Boost libs/numeric/odeint/examples/solar_system.cpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 solar system example for Hamiltonian stepper

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#include <bits/stdc++.h>
#include <bits/stdtr1c++.h>

#include <boost/circular_buffer.hpp>
#include <boost/numeric/odeint.hpp>

#include "solar_system.hpp"

#define tab "\t" 

using namespace std;
using namespace boost::numeric::odeint;


const size_t n = 3;
typedef point< double ,3 > point_type;
typedef std::tr1::array< point_type , n > state_type;
typedef std::tr1::array< double , n > mass_type;

typedef hamiltonian_stepper_rk< state_type > stepper_type;
typedef boost::circular_buffer< point_type > buffer_type;


const mass_type masses = {{ 1.0e9 , 1.0e9 , 1.0e12}};
const double gravitational_constant = 6.657e-20;

ostream& operator<<( ostream &out , const state_type &x )
{
    typedef state_type::value_type value_type;
    copy( x.begin() , x.end() ,
	  ostream_iterator< value_type >( out , "\n" ) );
    return out;
}

point_type get_mean( const state_type &x )
{
    point_type mean( 0.0 );
    if( x.empty() ) return mean;
    for( size_t i=0 ; i<x.size() ; ++i ) mean += x[i];
    mean /= double( x.size() );
    return mean;
}

point_type get_center_of_mass( const state_type &x ,
			       const mass_type &m )
{
    point_type mean( 0.0 );
    if( x.empty() ) return mean;
    double overall_mass = 0.0;
    for( size_t i=0 ; i<x.size() ; ++i )
    {
	overall_mass += m[i];
	mean += m[i] * x[i];
    }
    mean /= overall_mass;
    return mean;

}

void center_system( state_type &x , point_type mean )
{
    for( size_t i=0 ; i<x.size() ; ++i ) x[i] -= mean;
}


void solar_system( state_type &q , state_type &dpdt )
{
    point_type diff , tmp;
    fill( dpdt.begin() , dpdt.end() , 0.0 );
    for( size_t i=0 ; i<n ; ++i )
    {
	for( size_t j=i+1 ; j<n ; ++j )
	{
	    diff = q[j] - q[i];
	    tmp = gravitational_constant * diff / pow( abs( diff ) , 3.0 );
	    dpdt[i] += tmp * masses[j];
	    dpdt[j] -= tmp * masses[i];
	}
    }
}


int main( int argc , char **argv )
{
    state_type q , p;
    stepper_type stepper;

    fill( q.begin() , q.end() , 0.0 );
    fill( p.begin() , p.end() , 0.0 );
    q[0] = point_type( 0.0 , 1.0 , 0.0 );
    p[0] = point_type( 0.00001 , 0.0 , 0.0 );
    q[2] = point_type( 1.0 , 0.0 , 0.0 );

    center_system( q , get_center_of_mass( q , masses ) );
    center_system( p , get_center_of_mass( p , masses ) );

    const double dt = 1.0;

    buffer_type buffers[n];
    for( size_t i=0 ; i<n ; ++i ) buffers[i].set_capacity( 100 );

    cout << "unset key\n";
    const size_t ilen = 1000;
    double t = 0.0;
    while( true )
    {
	clog << get_mean( p ) << tab << get_mean( q ) << endl;
	for( size_t i=0 ; i<n ; ++i ) buffers[i].push_back( q[i] );

	cout << "p [-20:20][-20:20] '-' u 1:2 w l\n";
	for( size_t i=0 ; i<n ; ++i )
	{
	    copy( buffers[i].begin() , buffers[i].end() ,
		  ostream_iterator<point_type>( cout , "\n" ) );
	    cout << "\n";
	}
	cout << "e" << endl;

//	getchar();

//	cout << "p [-2:2][-2:2] '-' u 1:2 pt 7 ps 5 \n" << q << "e" << endl;

	for( size_t ii=0 ; ii<ilen ; ++ii,t+=dt )
	{
	    stepper.do_step( solar_system , q , p , dt );
	}
    }

    return 0;
}
