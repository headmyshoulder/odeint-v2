/*
 * integrate_functions.cpp
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#include <fstream>
#include <iostream>
#include <boost/array.hpp>
#include <vector>

#include <boost/numeric/odeint/stepper/implicit_euler.hpp>
#include <boost/numeric/odeint/stepper/rosenbrock4.hpp>
#include <boost/numeric/odeint/stepper/rosenbrock4_controller.hpp>

#include <boost/numeric/odeint/stepper/euler.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/dense_output_runge_kutta.hpp>

#include <boost/numeric/odeint/integrate/integrate.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>


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

typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::array< double , 3 > state_type;

std::ostream& operator<<( std::ostream &out , const vector_type &x )
{
	if( x.size() != 0 ) out << x[0];
	for( size_t i=1 ; i<x.size() ; ++i )
		out << " " << x[i];
	return out;
}

std::ostream& operator<<( std::ostream &out , const state_type &x )
{
	if( x.size() != 0 ) out << x[0];
	for( size_t i=1 ; i<x.size() ; ++i )
		out << " " << x[i];
	return out;
}




using namespace boost::numeric::odeint;

struct tmp_func
{
	std::ostream &m_out;
	tmp_func( std::ostream &out ) : m_out( out ) { }

	template< class T1 , class T2 >
	void operator()( const T1 &t1 , const T2 &t2 ) const
	{
		m_out << t2 << " " << t1 << "\n";
	}
};

int main( int argc , char **argv )
{
	using namespace std;

	state_type x1 = { { 10.0 , 10.0 , 10.0 } };
	vector_type x2( 3 , 10.0 );

//	integrate( implicit_euler< double >() , make_pair( lorenz() , lorenz_jacobi() ) , x2 , 0.0 , 10.0 , 0.1 , do_nothing_observer() );
//	integrate_n_steps( implicit_euler< double >() , make_pair( lorenz() , lorenz_jacobi() ) , x2 , 0.0 , 1000 , 0.1 );
//	integrate_adaptive( implicit_euler< double >() , make_pair( lorenz() , lorenz_jacobi() ) , x2 , 0.0 , 10.0 , 0.1 );

//	integrate( rosenbrock4< double >() , make_pair( lorenz() , lorenz_jacobi() ) , x2 , 0.0 , 10.0 , 0.1 );
//	integrate_n_steps( rosenbrock4< double >() , make_pair( lorenz() , lorenz_jacobi() ) , x2 , 0.0 , 1000 , 0.1 );
//	integrate_adaptive( rosenbrock4< double >() , make_pair( lorenz() , lorenz_jacobi() ) , x2 , 0.0 , 10.0 , 0.1 );
//
//	integrate( rosenbrock4_controller< rosenbrock4< double > >() , make_pair( lorenz() , lorenz_jacobi() ) , x2 , 0.0 , 10.0 , 0.1 );
//	integrate_n_steps( rosenbrock4_controller< rosenbrock4< double > >() , make_pair( lorenz() , lorenz_jacobi() ) , x2 , 0.0 , 1000 , 0.1 );
//	integrate_adaptive( rosenbrock4_controller< rosenbrock4< double > >() , make_pair( lorenz() , lorenz_jacobi() ) , x2 , 0.0 , 10.0 , 0.1 );


	integrate_const( euler< state_type >() , lorenz() , x1 , 0.0 , 0.1001 , 0.01 , tmp_func( cout ) );
//	integrate_n_steps( explicit_euler< state_type >() , lorenz() , x1 , 0.0 , 0.1 , 100 , cout << _1 << "\n" );
//	integrate_adaptive( explicit_euler< state_type >() , lorenz() , x1 , 0.0 , 10.0 , 0.1 , cout << _1 << "\n" );





	{
		// works
		ofstream fout( "integrate_controlled_rk54.dat" );
		size_t num_of_steps = integrate_const(
				controlled_runge_kutta< runge_kutta_cash_karp54< state_type > >() ,
				lorenz() , x1 , 0.0 , 50.0 , 0.1 , tmp_func( fout ) );
		clog << "Integrate controlled error stepper rk54 " << num_of_steps << endl;
	}



	{
		// seem to work, check
		typedef runge_kutta_dopri5< state_type > dopri5_type;
		typedef controlled_runge_kutta< dopri5_type > controlled_error_stepper_type;
		typedef dense_output_runge_kutta< controlled_error_stepper_type > stepper_type;

		controlled_error_stepper_type controlled_stepper( default_error_checker< double >( 1.0e-1 , 0.1 ) );

		ofstream fout( "integrate_controlled_dopri5.dat" );
		size_t num_of_steps = integrate_const(
				stepper_type( controlled_stepper ) , lorenz() , x1 , 0.0 , 50.0 , 0.001 , tmp_func( fout ) );
		clog << "Integrate denseoutput controlled dopri5 " << num_of_steps << endl;
	}

	{
		// seem to work, check
		typedef runge_kutta_dopri5< state_type > dopri5_type;
		typedef controlled_runge_kutta< dopri5_type > controlled_error_stepper_type;
		typedef dense_output_runge_kutta< controlled_error_stepper_type > stepper_type;

		controlled_error_stepper_type controlled_stepper( default_error_checker< double >( 1.0e-1 , 0.1 ) );

		ofstream fout( "integrate_adpative_controlled_dopri5.dat" );
		size_t num_of_steps = integrate_adaptive(
				stepper_type( controlled_stepper ) , lorenz() , x1 , 0.0 , 50.0 , 0.001 , tmp_func( fout ) );
		clog << "Integrate adaptive denseoutput controlled dopri5 " << num_of_steps << endl;
	}




	return true;
}
