/* Boost stepper_euler.cpp test file

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 This file tests the use of the euler stepper

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <vector>
#include <cmath>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/gsl_vector_adaptor.hpp>

using namespace boost::numeric::odeint;


void constant_system2( const gsl_vector &x , gsl_vector &dxdt , double t ) { gsl_vector_set( &dxdt , 0 , 1.0 ); }

const double eps = 1.0e-14;

int main( int argc , char **argv )
{

    gsl_vector *x_vec = gsl_vector_alloc( 1 );
    gsl_vector &x = *x_vec;

    explicit_euler< gsl_vector > euler;
    euler.do_step( constant_system2 , x , 0.0 , 0.1 );

    gsl_vector_free( x_vec );
}




