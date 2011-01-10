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
    : m_rb4() , m_atol( atol ) , m_rtol( rtol ) , m_first_step( true )
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
		static const time_type safe = 0.9 , fac1 = 5.0 , fac2 = 1.0 / 6.0;

		const size_t n = x.size();
		state_type xnew( n ) , xerr( n );
		m_rb4.do_step( sys , jacobi , x , t , xnew , dt , xerr );

		time_type err = error( xnew , xold , xerr );
		if ( err <= 1.0 )
		{
			if( m_first_step )
			{
				m_first_step = false;
			}
			else
			{
				time_type fac_pred = ( dt_old / dt ) * pow( err * err / errold , 0.25 ) / safe;
				fac_pred = std::max( fac2 , std::min( fac1 , fac_pred ) );
				fac = std::max( fac , fac_pred );
				dt_new = dt / fac;
			}

			dt_old = dt;
			errold = std::max( 0.01 , err );
			if( reject )
				dt_new = ( dt >= 0.0 ? std::min( dt_new , dt ) : std::max( dt_new , dt ) );
			dt_next = dt_new;
			reject=false;
		}
		else
		{
			time_type fac = std::max( fac2 ,std::min( fac1 , std::pow( err , 0.25 ) / safe ) );
			dt /= fac;
			return step_size_decreased;
		}

		return success_step_size_unchanged;
	}



private:

	rosenbrock4< time_type > m_rb4;
	time_type m_atol , m_rtol;
	bool m_first_step;
};


}
}
}


#endif /* ROSENBROCK4_HPP_ */
