/*
 * rosenbrock4.hpp
 *
 *  Created on: Jan 9, 2011
 *      Author: karsten
 */

#ifndef ROSENBROCK4_HPP_
#define ROSENBROCK4_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>


namespace boost {
namespace numeric {
namespace odeint {



template< class Value >
class rosenbrock4
{
public:

	typedef Value time_type;
    typedef boost::numeric::ublas::vector< time_type > state_type;
    typedef boost::numeric::ublas::matrix< time_type > matrix_type;
    typedef boost::numeric::ublas::permutation_matrix< size_t > pmatrix_type;

	rosenbrock4( void )
	{
	}


	template< class System , class Jacobi >
	void do_step( System &system , Jacobi &jacobi , const state_type &x , time_type t , state_type &xout , time_type dt , state_type &xerr )
	{
		const double gamma = 1.0;
		const double d1 = 1.0 , d2 = 1.0 , d3 = 1.0 , d4 = 1.0;
		const double c2 = 1.0 , c3 = 1.0 , c4 = 1.0;
		const double c21 = 1.0;
		const double a21 = 1.0;
		const double c31 = 1.0 , c32 = 1.0;
		const double a31 = 1.0 , a32 = 1.0;
		const double c41 = 1.0 , c42 = 1.0 , c43 = 1.0;
		const double a41 = 1.0 , a42 = 1.0 , a43 = 1.0;
		const double c51 = 1.0 , c52 = 1.0 , c53 = 1.0 , c54 = 1.0 ;
		const double a51 = 1.0 , a52 = 1.0 , a53 = 1.0 , a54 = 1.0 ;
		const double c61 = 1.0 , c62 = 1.0 , c63 = 1.0 , c64 = 1.0 , c65 = 1.0;


		const size_t n = x.size();
		matrix_type jac( n , n );
		pmatrix_type pm( n );
		state_type dfdt( n ) , dxdt( n );
		system( x , dxdt , t );
		jacobi( x , jac , t , dfdt );

		state_type g1( n ) , g2( n ) , g3( n ) , g4( n ) , g5( n );
		state_type xtmp( n ) , dxdtnew( n );


		jac *= -1.0;
		jac += 1.0 / gamma / dt * boost::numeric::ublas::identity_matrix< time_type >( n );
        boost::numeric::ublas::lu_factorize( jac , pm );

        for( size_t i=0 ; i<n ; ++i )
        	g1[i] = dxdt[i] + dt * d1 * dfdt[i];
        boost::numeric::ublas::lu_substitute( jac , pm , g1 );


        for( size_t i=0 ; i<n ; ++i )
        	xtmp[i] = x[i] + a21 * g1[i];
        system( xtmp , dxdtnew , t + c2 * dt );
        for( size_t i=0 ; i<n ; ++i )
        	g2[i] = dxdtnew[i] + dt * d2 * dfdt[i] + c21 * g1[i] / dt;
        boost::numeric::ublas::lu_substitute( jac , pm , g2 );


        for( size_t i=0 ; i<n ; ++i )
        	xtmp[i] = x[i] + a31 * g1[i] + a32 * g2[i];
        system( xtmp , dxdtnew , t + c3 * dt );
        for( size_t i=0 ; i<n ; ++i )
        	g3[i] = dxdtnew[i] + dt * d3 * dfdt[i] + ( c31 * g1[i] + c32 * g2[i] ) / dt;
        boost::numeric::ublas::lu_substitute( jac , pm , g3 );


        for( size_t i=0 ; i<n ; ++i )
        	xtmp[i] = x[i] + a41 * g1[i] + a42 * g2[i] + a43 * g3[i];
        system( xtmp , dxdtnew , t + c4 * dt );
        for( size_t i=0 ; i<n ; ++i )
        	g4[i] = dxdtnew[i] + dt * d4 * dfdt[i] + ( c41 * g1[i] + c42 * g2[i] + c43 * g3[i] ) / dt;
        boost::numeric::ublas::lu_substitute( jac , pm , g4 );


        for( size_t i=0 ; i<n ; ++i )
        	xtmp[i] = x[i] + a51 * g1[i] + a52 * g2[i] + a53 * g3[i] + a54 * g4[i];
        system( xtmp , dxdtnew , t + dt );
        for( size_t i=0 ; i<n ; ++i )
        	g5[i] = dxdtnew[i] + ( c51 * g1[i] + c52 * g2[i] + c53 * g3[i] + c54 * g4[i] ) / dt;
        boost::numeric::ublas::lu_substitute( jac , pm , g5 );

        for( size_t i=0 ; i<n ; ++i )
        	xtmp[i] += g5[i];
        system( xtmp , dxdtnew , t + dt );
        for( size_t i=0 ; i<n ; ++i )
        	xerr[i] = dxdtnew[i] + ( c61 * g1[i] + c62 * g2[i] + c63 * g3[i] + c64 * g4[i] + c65 * g5[i] ) / dt;
        boost::numeric::ublas::lu_substitute( jac , pm , xerr );

        for( size_t i=0 ; i<n ; ++i )
        	xout[i] = xtmp[i] + xerr[i];
	}

	template< class System , class Jacobi >
	void do_step( System &system , Jacobi &jacobi , state_type &x , time_type t , time_type dt , state_type &xerr )
	{
		state_type out( x.size() );
		do_step( system , jacobi , x , t , out , dt , xerr );
		x = out;
	}


private:

};



template< class Value >
class rosenbrock4_controller
{
public:

	typedef Value time_type;
    typedef boost::numeric::ublas::vector< time_type > state_type;
    typedef boost::numeric::ublas::matrix< time_type > matrix_type;
    typedef boost::numeric::ublas::permutation_matrix< size_t > pmatrix_type;

	rosenbrock4_controller( time_type atol = 1.0e-6 , time_type rtol = 1.0e-6 )
    : m_rb4() , m_atol( atol ) , m_rtol( rtol )
	{
	}

	time_type error( const state_type &x , const state_type &xold , const state_type &xerr )
	{
		const size_t n = x.size();
		time_type err = 0.0 , sk = 0.0;
		for( size_t i=0 ; i<n ; ++i )
		{
			sk = m_atol + m_rtol * std::max( std::abs( x[i] ) , std::abs( xold[i] ) );
			err += xerr[i] * xerr[i] / sk / sk;
		}
		return sqrt( err / time_type( n ) );
	}


	template< class System , class Jacobi >
	boost::numeric::odeint::controlled_step_result
	try_step( System &sys , Jacobi &jacobi , state_type &x , time_type &t , time_type &dt )
	{
		// hier gehts weiter
//		const size_t n = x.size();
		return success_step_size_unchanged;
	}

//	template <class D> stepperross.h
//	StepperRoss<D>::Controller::Controller() : reject(false), first_step(true) {}
//	Step size controller for fourth-order Rosenbrock method.
//	template <class D>
//	bool StepperRoss<D>::Controller::success(Doub err, Doub &h) {
//	Returns true if err  1, false otherwise. If step was successful, sets hnext to the estimated
//	optimal stepsize for the next step. If the step failed, reduces h appropriately for another try.
//	static const Doub safe=0.9,fac1=5.0,fac2=1.0/6.0;
//	Doub fac=MAX(fac2,MIN(fac1,pow(err,0.25)/safe));
//	Doub hnew=h/fac; Ensure 1=fac1  hnew=h  1=fac2.
//	if (err <= 1.0) { Step succeeded.
//	if (!first_step) { Predictive control.
//	Doub facpred=(hold/h)*pow(err*err/errold,0.25)/safe;
//	facpred=MAX(fac2,MIN(fac1,facpred));
//	fac=MAX(fac,facpred);
//	hnew=h/fac;
//	}
//	first_step=false;
//	hold=h;
//	errold=MAX(0.01,err);
//	if (reject) Don’t let step increase if last one was rejected.
//	hnew=(h >= 0.0 ? MIN(hnew,h) : MAX(hnew,h));
//	hnext=hnew;
//	reject=false;
//	return true;
//	} else { Truncation error too large, reduce stepsize.
//	h=hnew;
//	reject=true;
//	return false;
//	}
//	}


private:

	rosenbrock4< time_type > m_rb4;
	time_type m_atol , m_rtol;

};


}
}
}


#endif /* ROSENBROCK4_HPP_ */
