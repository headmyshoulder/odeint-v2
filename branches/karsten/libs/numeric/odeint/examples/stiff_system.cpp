/*
 * rosenbrock4.cpp
 *
 *  Created on: Jan 9, 2011
 *      Author: karsten
 */

#include <iostream>
#include <fstream>
#include <utility>
#include <tr1/array>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

//[ stiff_system_definition
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

struct stiff_system
{
	void operator()( const vector_type &x , vector_type &dxdt , double t )
	{
		dxdt[ 0 ] = -101.0 * x[ 0 ] - 100.0 * x[ 1 ];
		dxdt[ 1 ] = x[ 0 ];
	}
};

struct stiff_system_jacobi
{
	void operator()( const vector_type &x , matrix_type &J , const double &t , vector_type &dfdt )
	{
		J( 0 , 0 ) = -101.0;
		J( 0 , 1 ) = -100.0;
		J( 1 , 0 ) = 1.0;
		J( 1 , 1 ) = 0.0;
		dfdt[0] = 0.0;
		dfdt[1] = 0.0;
	}
};
//]


/*
//[ stiff_system_alternative_definition
struct stiff_system
{
	template< class State >
	void operator()( const State &x , State &dxdt , double t )
	{
		...
	}
};

struct stiff_system_jacobi
{
	template< class State , class Matrix >
	void operator()( const State &x , Matrix &J , const double &t , State &dfdt )
	{
		...
	}
};
//]
*/





void rk54_ck_controlled_with_stiff_system( void )
{
	typedef explicit_error_rk54_ck< vector_type > stepper_type;
	typedef controlled_error_stepper< stepper_type > controlled_stepper_type;
	controlled_stepper_type stepper;

	vector_type x( 3 , 1.0 );
	double t = 0.0 , dt = 0.00001;
	ofstream fout( "rk54_ck_controller_stiff.dat" );
	size_t count = 0;
	while( t < 50.0 )
	{
		fout << t << "\t" << dt << "\t" << stepper.last_error() << "\t";
		fout << x[0] << "\t" << x[1] << "\t" << x[2] << "\t";
		fout <<std::endl;

		size_t trials = 0;
		while( trials < 100 )
		{
			if( stepper.try_step( stiff_system() , x , t , dt ) !=  step_size_decreased )
				break;
			++trials;
		}
		if( trials == 100 )
		{
			cerr << "Error : stepper did not converge! " << endl;
			break;
		}
		++count;
	}
	clog << "RK 54 : " << count << endl;
}






int main( int argc , char **argv )
{
//[ integrate_stiff_system
	typedef rosenbrock4< double > stepper_type;
	typedef rosenbrock4_controller< stepper_type > controlled_stepper_type;

	vector_type x( 3 , 1.0 ) , xerr( 3 );

	size_t num_of_steps = integrate_const( controlled_stepper_type() ,
			make_pair( stiff_system() , stiff_system_jacobi() ) ,
			x , 0.0 , 50.0 , 0.01 );
//]
	clog << num_of_steps << endl;

	rk54_ck_controlled_with_stiff_system();

	return 0;
}
