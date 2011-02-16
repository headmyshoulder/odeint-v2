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

#include <boost/numeric/odeint/stepper/rosenbrock4.hpp>
#include <boost/numeric/odeint/stepper/rosenbrock4_controller.hpp>
#include <boost/numeric/odeint/stepper/explicit_error_rk54_ck.hpp>
#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;


struct lorenz
{
	template< class StateType >
	void operator()( const StateType &x , StateType &dxdt , const double &t )
	{

		dxdt[0] = sigma * ( x[1] - x[0] );
		dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
		dxdt[2] = x[0] * x[1] - b * x[2];
	}
};

struct lorenz_jacobi
{
	template< class State , class Matrix >
	void operator()( const State &x , Matrix &J , const double &t , State &dfdt )
	{
		J( 0 , 0 ) = -sigma;
		J( 0 , 1 ) = sigma;
		J( 0 , 2 ) = 0.0;
		J( 1 , 0 ) = R - x[2];
		J( 1 , 1 ) = -1.0;
		J( 1 , 2 ) = -x[0];
		J( 2 , 0 ) = x[1];
		J( 2 , 1 ) = x[0];
		J( 2 , 2 ) = -b;

		dfdt[0] = 0.0;
		dfdt[1] = 0.0;
		dfdt[2] = 0.0;
	}
};


struct stiff_system
{
	template< class State >
	void operator()( const State &x , State &dxdt , double t )
	{
		dxdt[ 0 ] = -101.0 * x[ 0 ] - 100.0 * x[ 1 ];
		dxdt[ 1 ] = x[ 0 ];
	}
};

struct stiff_system_jacobi
{
	template< class State , class Matrix >
	void operator()( const State &x , Matrix &J , const double &t , State &dfdt )
	{
		J( 0 , 0 ) = -101.0;
		J( 0 , 1 ) = -100.0;
		J( 1 , 0 ) = 1.0;
		J( 1 , 1 ) = 0.0;
		dfdt[0] = 0.0;
		dfdt[1] = 0.0;
	}
};





void test_rosenbrock_stepper_with_lorenz( void )
{
	const static size_t dim = 3;
	typedef rosenbrock4< double > stepper_type;
	typedef stepper_type::state_type state_type;
	typedef stepper_type::matrix_type matrix_type;

	const double dt = 0.01;
	size_t steps = 1000;
	double x0 = -12.0 , y0 = -12.0 , z0 = 20.0;
	state_type x( dim ) , xerr( dim );
	double t = 0.0;

	stepper_type stepper;
	x[0] = x0 ; x[1] = y0 ; x[2] = z0;

	ofstream fout( "rosenbrock_stepper_lorenz.dat" );
	fout.precision( 14 );
	size_t count = 0;
	while( count < steps )
	{
		fout << t << " ";
		fout << x[0] << " " << x[1] << " " << x[2] << " ";
		fout << xerr[0] << " " << xerr[1] << " " << xerr[2] << " ";
		fout <<std::endl;

		stepper.do_step( make_pair( lorenz() , lorenz_jacobi() ) , x , t , dt , xerr );
		++count;
		t += dt;
	}
}

void rk54_ck_stepper_with_lorenz( void )
{
	typedef std::tr1::array< double , 3 > state_type2;
	typedef explicit_error_rk54_ck< state_type2 > stepper_type2;
	stepper_type2 rk_stepper;

	double x0 = -12.0 , y0 = -12.0 , z0 = 20.0;
	state_type2 x = {{ x0 , y0 , z0 }} , xerr = {{ 0.0 , 0.0 , 0.0 }};
	double t = 0.0;

	ofstream fout( "rk54_ck_stepper_with_lorenz.dat" );
	fout.precision( 14 );
	size_t count = 0;
	size_t steps = 1000;
	const double dt = 0.01;
	while( count < steps )
	{
		fout << t << "\t";
		fout << x[0] << "\t" << x[1] << "\t" << x[2] << "\t";
		fout << xerr[0] << "\t" << xerr[1] << "\t" << xerr[2] << "\t";
		fout <<std::endl;

		rk_stepper.do_step( lorenz() , x , t , dt , xerr );
		++count;
		t += dt;
	}
}

void test_controlled_rosenbrock_with_stiff_system( void )
{
	const static size_t dim = 3;
	typedef rosenbrock4< double > stepper_type;
	typedef stepper_type::state_type state_type;
	typedef stepper_type::matrix_type matrix_type;
	typedef rosenbrock4_controller< stepper_type > controlled_stepper_type;

	state_type x( dim ) , xerr( dim );
	double t = 0.0 , dt = 0.00001;

	controlled_stepper_type controlled_stepper;

	x[0] = 1.0 ; x[1] = 1.0 ; x[2] = 0.0;

	ofstream fout( "rosenbrock_controller_stiff.dat" );
	size_t count = 0;
	while( t < 50.0 )
	{
		fout << t << "\t" << dt << "\t" << controlled_stepper.last_error() << "\t";
		fout << x[0] << "\t" << x[1] << "\t" << x[2] << "\t";
		fout <<std::endl;

		size_t trials = 0;
		while( trials < 100 )
		{
			if( controlled_stepper.try_step( make_pair( stiff_system() , stiff_system_jacobi() ) , x , t , dt ) !=  step_size_decreased )
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
	clog << "Rosenbrock: " << count << endl;
}

void rk54_ck_controlled_with_stiff_system( void )
{
	typedef std::tr1::array< double , 3 > state_type2;
	typedef explicit_error_rk54_ck< state_type2 > stepper_type2;
	typedef controlled_error_stepper< stepper_type2 > controlled_stepper_type2;
	controlled_stepper_type2 stepper;

	state_type2 x = {{ 1.0 , 1.0 , 0.0 }};
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
	test_rosenbrock_stepper_with_lorenz();
	rk54_ck_stepper_with_lorenz();
	test_controlled_rosenbrock_with_stiff_system();
	rk54_ck_controlled_with_stiff_system();

	return 0;
}
