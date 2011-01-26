/* Boost implicit_euler.cpp example file

 Copyright 2010 Karsten Ahnert
 Copyright 2010 Mario Mulansky

 This file shows the use of the implicit euler stepper

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <iostream>
#include <utility>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/implicit_euler.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#define tab '\t'

typedef double value_type;
typedef boost::numeric::ublas::vector< value_type > state_type;
typedef boost::numeric::ublas::matrix< value_type > matrix_type;

void stiff_system( state_type &x , state_type &dxdt , value_type t )
{
	dxdt( 0 ) = -101.0 * x( 0 ) - 100.0 * x( 1 );
	dxdt( 1 ) = x( 0 );
}

void jacobi( state_type &x , matrix_type &jacobi , value_type t )
{
	jacobi( 0 , 0 ) = -101.0;
	jacobi( 0 , 1 ) = -100.0;
	jacobi( 1 , 0 ) = 1.0;
	jacobi( 1 , 1 ) = 0.0;
}

using namespace std;
using namespace boost::numeric::odeint;

int main( void )
{
	explicit_euler< state_type , value_type /*, vector_space_algebra */> expl_euler;
	implicit_euler< value_type > impl_euler;

	state_type x1( 2 );
	x1( 0 ) = 1.0; x1( 1 ) = 0.0;
	state_type x2( x1 );

	const value_type dt = 0.01; // for dt >= 0.01 the euler method gets unstable
	const size_t steps = 1000;

	value_type t = 0.0;
	for( size_t step = 0 ; step < steps ; ++step , t+=dt )
	{
		clog << step << " of " << steps << endl;
		cout << t << tab << x1( 0 ) << tab << x1( 1 ) << tab << x2( 0 ) << tab << x2( 1 ) << tab;
		cout << 1.0 / 99.0 * exp( -100.0 * t ) * ( 100.0  - exp( 99.0 * t ) ) << tab;
		cout << 1.0 / 99.0 * exp( -100.0 * t ) * ( -1.0 + exp( 99.0 * t ) ) << endl;

		expl_euler.do_step( stiff_system , x1 , t , dt );
		impl_euler.do_step( make_pair( stiff_system , jacobi ) , x2 , t , dt );
	}

}
