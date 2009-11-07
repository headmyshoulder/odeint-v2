/* Boost numeric/odeint/examples/lorenz_stepper_cmp.cpp
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Shows the usage of odeint, and integrates the Lorenz equations,
 a deterministic chaotic system

 dx/dt = sigma * ( x - y)
 dy/dt = R*x - y - x*z
 dz/dt = x*y - b*z

 with sigma = 10, r=28, b = 8/3

 Furhtermore, it compares the solution with the gsl rk4 routine.

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <tr1/array>

#include <boost/numeric/odeint.hpp>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>


#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;


typedef std::tr1::array< double , 3 > state_type;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

void lorenz( state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

int lorenz_gsl( double t , const double y[] , double f[] , void *params)
{
    f[0] = sigma * ( y[1] - y[0] );
    f[1] = R * y[0] - y[1] - y[0] * y[2];
    f[2] = y[0]*y[1] - b * y[2];
    return GSL_SUCCESS;
}

void lorenz2( double *x , double *dxdt )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

void rk4_lorenz( double *x , double h )
{
    const size_t n = 3;

    double hh , h6 , dxdt[n] , dxm[n] , dxt[n] , xt[n];
    hh = h*0.5;
    h6 = h/6.0;

    lorenz2( x , dxdt );
    for( size_t i=0;i<n;i++) xt[i]=x[i]+hh*dxdt[i];
    lorenz2(xt,dxt);
    for ( size_t i=0;i<n;i++) xt[i]=x[i]+hh*dxt[i];
    lorenz2(xt,dxm);
    for( size_t i=0;i<n;i++)
    {
        xt[i]=x[i]+h*dxm[i];
        dxm[i] += dxt[i];
    }
    lorenz2(xt,dxt);
    for ( size_t i=0;i<n;i++)
        x[i] += h6*(dxdt[i]+dxt[i]+2.0*dxm[i]);
}


int main( int argc , char **argv )
{
    const double dt = 0.01;
    const size_t olen = 10000000;


    state_type x_start = {{ -13.90595 , -15.54127 , 32.48918 }};

    clock_t start , end;
    double t;



    // odeint method
    state_type x1 = x_start , x1_err;
    ode_step_half_step< ode_step_runge_kutta_4< state_type > > stepper;

    start= clock();
    t = 0.0;
    for( size_t oi=0 ; oi<olen ; ++oi,t+=dt )
	stepper.next_step( lorenz , x1 , t , dt , x1_err );
    end = clock();
    clog << "odeint array : " << double ( end - start ) / double( CLOCKS_PER_SEC ) << endl;




    // gsl method
    gsl_odeiv_step *s = gsl_odeiv_step_alloc( gsl_odeiv_step_rk4 , 3);
    gsl_odeiv_system sys = { lorenz_gsl , 0 , 3 , 0 };
    double x2[3] = { x_start[0] , x_start[1] , x_start[2] } , x2_err[3];

    start= clock();
    t = 0.0;
    for( size_t oi=0 ; oi<olen ; ++oi,t+=dt )
	gsl_odeiv_step_apply ( s , t , dt , x2 , x2_err , 0 , 0 , &sys );
    end = clock();
    clog << "gsl rk4 : " << double ( end - start ) / double( CLOCKS_PER_SEC ) << endl;
     


    // just check, if rk4 is working correctly
    const size_t tslen = 10000;
    x1 = x_start;
    copy( x_start.begin() , x_start.end() , x2 );
    double x3[3] = { x_start[0] , x_start[1] , x_start[2] };
    cout.precision( 14 );
    cout.flags( ios::scientific );

    t = 0.0;
    for( size_t i=0 ; i<tslen ; ++i,t+=dt )
    {
	stepper.next_step( lorenz , x1 , t , dt );
	gsl_odeiv_step_apply ( s , t , dt , x2 , x2_err , 0 , 0 , &sys );
	rk4_lorenz( x3 , dt );
	cout << t << tab << x1[0] << tab << x2[0] << tab << x3[0] << endl;
    }

    gsl_odeiv_step_free (s);


    return 0;
}



/*
  Compile with
  g++ -Wall -O3 -I$BOOST_ROOT -I../../../../ lorenz_stepper_cmp.cpp -lgsl -lgslcblas
*/
