/*
 * integrate_functions.cpp
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#include <tr1/array>

#include <boost/numeric/odeint/stepper/implicit_euler.hpp>
#include <boost/numeric/odeint/stepper/rosenbrock4.hpp>
#include <boost/numeric/odeint/stepper/rosenbrock4_controller.hpp>

#include <boost/numeric/odeint/stepper/explicit_euler.hpp>
#include <boost/numeric/odeint/stepper/explicit_error_dopri5.hpp>
#include <boost/numeric/odeint/stepper/explicit_error_rk54_ck.hpp>
#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>
#include <boost/numeric/odeint/stepper/dense_output_explicit.hpp>
#include <boost/numeric/odeint/stepper/dense_output_controlled_explicit_fsal.hpp>

#include <boost/numeric/odeint/integrate/integrate.hpp>



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
typedef std::tr1::array< double , 3 > state_type;

using namespace std;
using namespace boost::numeric::odeint;

int main( int argc , char **argv )
{
	state_type x1;
	vector_type x2( 3 );

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


//#include <boost/numeric/odeint/stepper/explicit_euler.hpp>
//#include <boost/numeric/odeint/stepper/explicit_error_dopri5.hpp>
//#include <boost/numeric/odeint/stepper/explicit_error_rk54_ck.hpp>
//#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>
//#include <boost/numeric/odeint/stepper/dense_output_explicit.hpp>
//#include <boost/numeric/odeint/stepper/dense_output_controlled_explicit_fsal.hpp>


	return true;
}
