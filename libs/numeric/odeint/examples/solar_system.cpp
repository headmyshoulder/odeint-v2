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

#include <boost/numeric/odeint/container_traits_tr1_array.hpp>
#include <boost/numeric/odeint/hamiltonian_stepper_euler.hpp>
#include "point_type.hpp"

#define tab "\t" 

using namespace std;
using namespace boost::numeric::odeint;

// we simulate 5 planets and the sun
const size_t n = 6;

//[ state_type_definition
typedef point< double , 3 > point_type;
typedef std::tr1::array< point_type , n > container_type;
typedef std::tr1::array< double , n > mass_type;
typedef pair< container_type , container_type > state_type;
//]



const double gravitational_constant = 2.95912208286e-4;


ostream& operator<<( ostream &out , const container_type &x )
{
    typedef container_type::value_type value_type;
    copy( x.begin() , x.end() ,
	  ostream_iterator< value_type >( out , "\n" ) );
    return out;
}


point_type center_of_mass( const container_type &x , const mass_type &m )
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

point_type center( const container_type &x )
{
    point_type mean( 0.0 );
    for( size_t i=0 ; i<x.size() ; ++i ) mean += x[i];
    if( !x.empty() ) mean /= x.size();
    return mean;
}


double energy( const container_type &q , const container_type &p , const mass_type &masses )
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


//[ coordinate_function
struct solar_system
{
    const mass_type &m_masses;

    solar_system( const mass_type &masses ) : m_masses( masses ) { }

    void operator()( const container_type &q , container_type &dpdt ) const
    {
        const size_t n = q.size();
        fill( dpdt.begin() , dpdt.end() , point_type( 0.0 , 0.0 , 0.0 ) );
        for( size_t i=0 ; i<n ; ++i )
        {
            for( size_t j=i+1 ; j<n ; ++j )
            {
                point_type diff = q[j] - q[i];
                double d = abs( diff );
                diff *= ( gravitational_constant / d / d / d * m_masses[i] * m_masses[j] ) ;
                dpdt[i] += diff ;
                dpdt[j] -= diff ;
            }
        }
    }
};
//]

//[ momentum_function
struct solar_system_momentum
{
    const mass_type &m_masses;

    solar_system_momentum( const mass_type &masses ) : m_masses( masses ) { }

    void operator()( const container_type &p , container_type &dqdt ) const
    {
        for( size_t i=0 ; i<p.size() ; ++i ) dqdt[i] = p[i] / m_masses[i];
    }
};
//]



int main( int argc , char **argv )
{
    const mass_type masses = {{
            1.00000597682 ,      // sun
            0.000954786104043 ,  // jupiter
            0.000285583733151 ,  // saturn
            0.0000437273164546 , // uranus
            0.0000517759138449 , // neptune
            1.0 / ( 1.3e8 )      // pluto
        }};

    container_type q = {{
            point_type( 0.0 , 0.0 , 0.0 ) ,                        // sun
            point_type( -3.5023653 , -3.8169847 , -1.5507963 ) ,   // jupiter
            point_type( 9.0755314 , -3.0458353 , -1.6483708 ) ,    // saturn
            point_type( 8.3101420 , -16.2901086 , -7.2521278 ) ,   // uranus
            point_type( 11.4707666 , -25.7294829 , -10.8169456 ) , // neptune
            point_type( -15.5387357 , -25.2225594 , -3.1902382 )   // pluto
        }};

    container_type p = {{
            point_type( 0.0 , 0.0 , 0.0 ) ,                        // sun
            point_type( 0.00565429 , -0.00412490 , -0.00190589 ) , // jupiter
            point_type( 0.00168318 , 0.00483525 , 0.00192462 ) ,   // saturn
            point_type( 0.00354178 , 0.00137102 , 0.00055029 ) ,   // uranus
            point_type( 0.00288930 , 0.00114527 , 0.00039677 ) ,   // neptune
            point_type( 0.00276725 , -0.00170702 , -0.00136504 )   // pluto
        }};

    // remove center of mass velocity
    point_type com = center_of_mass( p , masses );
    for( size_t i=0 ; i<n ; ++i ) p[i] -= com;
    for( size_t i=0 ; i<n ; ++i ) p[i] *= masses[i];


    const double dt = 1.0;
    double t = 0.0;

//[ integration_solar_system
    typedef hamiltonian_stepper_euler< container_type > stepper_type;
    stepper_type stepper;
    state_type state = make_pair( q , p );
    for( size_t c = 0 ; c<2000 ; ++c )
    {
        clog << t << tab << energy( state.first , state.second , masses ) << tab;
        clog << center_of_mass( state.first , masses ) << tab;
        clog << center_of_mass( state.second , masses ) << endl;

        cout << t;
        for( size_t i=0 ; i<n ; ++i ) cout << tab << state.first[i];
        cout << endl;

        for( size_t i=0 ; i<100 ; ++i,t+=dt )
            stepper.do_step( make_pair( solar_system_momentum( masses ) ,
                                        solar_system( masses ) ) ,
                             state , t , dt );
        t += dt;
    }
//]

    return 0;
}


/*
Plot with gnuplot:
p "solar_system.dat" u 2:4 w l,"solar_system.dat" u 5:7 w l,"solar_system.dat" u 8:10 w l,"solar_system.dat" u 11:13 w l,"solar_system.dat" u 14:16 w l,"solar_system.dat" u 17:19 w l
*/
