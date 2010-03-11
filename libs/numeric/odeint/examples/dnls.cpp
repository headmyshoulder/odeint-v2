/* Boost numeric/odeint/examples/lorenz_integrator.cpp
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 compares the steppers rk4 and rk78
 the system is the dnls, which is complex and Hamiltonian

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <iostream>
#include <vector>
#include <iterator>
#include <list>
#include <algorithm>
#include <tr1/array>

#include <boost/numeric/odeint.hpp>

#define tab "\t"

using namespace std;
using namespace std::tr1;
using namespace boost::numeric::odeint;


const size_t n = 64;
const double beta = 1.0;

typedef array< complex<double> , n > state_type;

const complex<double> II( 0.0 , -1.0 );

void dnls( state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = II * ( beta * norm( x[0] ) * x[0] + x[1] + x[n-1] );
    for( size_t i=1 ; i<n-1 ; ++i )
        dxdt[i] = II * ( beta * norm( x[i] ) * x[i] + x[i+1] + x[i-1] );
    dxdt[n-1] = II * ( beta * norm( x[n-1] ) * x[n-1] + x[0] + x[n-2] );
}

double norm( const state_type &x )
{
    double nn = 0.0;
    state_type::const_iterator iter = x.begin() ;
    while( iter != x.end() ) nn += norm( *iter++ );
    return nn;
}

double energy( const state_type &x )
{
    double e = 0.0 , nn;
    for( size_t i=0 ; i<n-1 ; ++i )
    {
        nn = norm( x[i] );
        e += 0.5 * beta * nn * nn + 2.0 * ( x[i]*conj(x[i+1]) ).real();
    }
    nn = norm( x[n-1] );
    e += 0.5 * beta * nn * nn + 2.0 * ( x[n-1]*conj(x[0]) ).real();
    return e;
}

ostream& operator<<( ostream &out , const state_type &x )
{
    state_type::const_iterator iter = x.begin() ;
    while( iter != x.end() )
    {
        const complex<double> &y = *iter++;
        out << y.real() << tab << y.imag() << tab << norm( y ) << "\n";
    }
    return out;
}

int main( int argc , char **argv )
{
    state_type x;

    generate( x.begin() , x.end() , drand48 );

    state_type x1( x ) , x2( x );
    stepper_rk4< state_type > rk4;
    stepper_rk78_fehlberg< state_type > rk78;

    double norm0 = norm( x1 ) , energy0 = energy( x1 );

    const size_t olen = 10000 , ilen = 100;
    const double dt = 0.01;
    double t = 0.0;
    cout.precision(14);
    cout.flags( ios::scientific );
    for( size_t oi=0 ; oi<olen ; ++oi )
    {
        double norm1 = norm( x1 ) , norm2 = norm( x2 );
        double energy1 = energy( x1 ) , energy2 = energy( x2 );

        cout << t << tab;
        cout << norm1 << tab << norm2 << tab;
        cout << energy1 << tab << energy2 << tab;
        cout << norm1 - norm0 << tab << norm2 - norm0 << tab;
        cout << energy1 - energy0 << tab << energy2 - energy0 << endl;

        for( size_t ii=0 ; ii<ilen ; ++ii,t+=dt )
        {
            rk4.do_step( dnls , x1 , t , dt );
            rk78.do_step( dnls , x2 , t , dt );
        }
    }



    return 0;
}
