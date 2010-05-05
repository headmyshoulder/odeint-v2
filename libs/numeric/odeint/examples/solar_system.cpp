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
#include <boost/ref.hpp>
#include <boost/numeric/odeint.hpp>

#include "point_type.hpp"

#define tab "\t" 

using namespace std;
using namespace boost::numeric::odeint;

// we simulate 5 planets and the sun
const size_t n = 6;

typedef point< double , 3 > point_type;
typedef std::tr1::array< point_type , n > state_type;
typedef std::tr1::array< double , n > mass_type;
typedef hamiltonian_stepper_euler< state_type > stepper_type;


const double gravitational_constant = 2.95912208286e-4;


ostream& operator<<( ostream &out , const state_type &x )
{
    typedef state_type::value_type value_type;
    copy( x.begin() , x.end() ,
	  ostream_iterator< value_type >( out , "\n" ) );
    return out;
}


point_type center_of_mass( const state_type &x , const mass_type &m )
{
    double overall_mass = 0.0;
    point_type mean( 0.0 );
    for( size_t i=0 ; i<x.size() ; ++i )
    {
	overall_mass += m[i];
	mean += m[i] * x[i];
    }
    if( x.size() != 0 ) mean /= overall_mass;
    return mean;
}


double energy( const state_type &q , const state_type &p ,
               const mass_type &masses )
{
    const size_t n = q.size();
    double en = 0.0;
    for( size_t i=0 ; i<n ; ++i )
    {
        en += 0.5 * norm( p[i] ) / masses[i];
	for( size_t j=0 ; j<i ; ++j )
        {
            double diff = abs( q[i] - q[j] );
            en -= gravitational_constant * masses[j] * masses[i] / diff;
        }
    }
    return en;
}

struct solar_system
{
    mass_type &m_masses;

    solar_system( mass_type &masses ) : m_masses( masses ) { }

    void operator()( const state_type &q , state_type &dpdt )
    {
        const size_t n = q.size();
        fill( dpdt.begin() , dpdt.end() , 0.0 );
        for( size_t i=0 ; i<n ; ++i )
        {
            for( size_t j=i+1 ; j<n ; ++j )
            {
                point_type diff = q[j] - q[i];
		double d = abs( diff );
                diff = gravitational_constant * diff / d / d / d;
                dpdt[i] += diff * m_masses[j];
                dpdt[j] -= diff * m_masses[i];
            }
        }
    }
};


int main( int argc , char **argv )
{
    mass_type masses = {{
            1.00000597682 ,      // sun
            0.000954786104043 ,  // jupiter
            0.000285583733151 ,  // saturn
            0.0000437273164546 , // uranus
            0.0000517759138449 , // neptune
            1.0 / ( 1.3e8 )      // pluto
        }};

    state_type q = {{
            point_type( 0.0 , 0.0 , 0.0 ) ,                        // sun
            point_type( -3.5023653 , -3.8169847 , -1.5507963 ) ,   // jupiter
            point_type( 9.0755314 , -3.0458353 , -1.6483708 ) ,    // saturn
            point_type( 8.3101420 , -16.2901086 , -7.2521278 ) ,   // uranus
            point_type( 11.4707666 , -25.7294829 , -10.8169456 ) , // neptune
            point_type( -15.5387357 , -25.2225594 , -3.1902382 )   // pluto
        }};

    state_type p = {{
            point_type( 0.0 , 0.0 , 0.0 ) ,                        // sun
            point_type( 0.00565429 , -0.00412490 , -0.00190589 ) , // jupiter
            point_type( 0.00168318 , 0.00483525 , 0.00192462 ) ,   // saturn
            point_type( 0.00354178 , 0.00137102 , 0.00055029 ) ,   // uranus
            point_type( 0.00288930 , 0.00114527 , 0.00039677 ) ,   // neptune
            point_type( 0.00276725 , -0.00170702 , -0.00136504 )   // pluto
        }};


    point_type qmean = center_of_mass( q , masses );
    point_type pmean = center_of_mass( p , masses );
    for( size_t i=0 ; i<n ; ++i ) { q[i] -= qmean ; p[i] -= pmean; }

    stepper_type stepper;

    const double dt = 100.0;
    double t = 0.0;
    while( t < 10000000.0 )
    {
	clog << t << tab << energy( q , p , masses ) << tab;
        clog << center_of_mass( q , masses ) << tab;
        clog << center_of_mass( p , masses ) << endl;

	cout << t;
        for( size_t i=0 ; i<n ; ++i ) cout << tab << q[i];
	cout << endl;

        for( size_t i=0 ; i<1 ; ++i,t+=dt )
	    stepper.do_step( solar_system( masses ) , q , p , dt );
        t += dt;
    }

    return 0;
}


/*
Plot with gnuplot:
p "solar_system.dat" u 2:4 w l,"solar_system.dat" u 5:7 w l,"solar_system.dat" u 8:10 w l,"solar_system.dat" u 11:13 w l,"solar_system.dat" u 14:16 w l,"solar_system.dat" u 17:19 w l
*/
