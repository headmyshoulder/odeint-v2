/*
 * bulirsch_stoer.cpp
 *
 *  Created on: Aug 8, 2011
 *      Author: mario
 */

#include <iostream>
#include <fstream>
#include <cmath>

#include <boost/array.hpp>

#include <boost/numeric/odeint/config.hpp>

#include <boost/lambda/lambda.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>

using namespace std;
using namespace boost::lambda;
using namespace boost::numeric::odeint;

typedef boost::array< double , 1 > state_type;

/*
 * x' = ( - x*sin t  + 2 tan x ) y
 * with x( pi/6 ) = 2/sqrt(3) the analytic solution is 1/cos t
 */

void rhs( const state_type &x , state_type &dxdt , const double t )
{
    dxdt[0] = ( - x[0] * sin( t ) + 2.0 * tan( t ) ) * x[0];
}

void rhs2( const state_type &x , state_type &dxdt , const double t )
{
    dxdt[0] = sin(t);
}


ofstream out;

void write_out( const state_type &x , const double t )
{
    out << t << '\t' << x[0] << endl;
}

int main()
{
    bulirsch_stoer_dense_out< state_type > stepper( 1E-10 , 0.0 , 0.0 , 0.0 );

    state_type x = {{ 0.0 }};

    //double t = M_PI/6.0;
    double t = 0.0;
    double dt = 0.01;
    //double t_end = M_PI/2.0 - 0.1;
    double t_end = 100.0;

    out.open( "bs.dat" );
    integrate_const( stepper , rhs2 , x , t , t_end , dt , write_out );
    out.close();

    x[0] = 0.0;

    out.open( "bs2.dat" );
    integrate_adaptive( stepper , rhs2 , x , t , t_end , dt , write_out );
    out.close();

    typedef runge_kutta_dopri5< state_type > dopri5_type;
    typedef controlled_error_stepper< dopri5_type > controlled_dopri5_type;
    typedef dense_output_controlled_explicit_fsal< controlled_dopri5_type > dense_output_dopri5_type;

    dense_output_dopri5_type dopri5( controlled_dopri5_type( dopri5_type() , default_error_checker< double >( 1E-10 , 0.0 , 0.0 , 0.0 )  ) );

    x[0] = 0.0;

    out.open( "bs3.dat" );
    integrate_adaptive( dopri5 , rhs2 , x , t , t_end , dt , write_out );
    out.close();

}
