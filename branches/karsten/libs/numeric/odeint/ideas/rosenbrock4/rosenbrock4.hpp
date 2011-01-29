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

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>


namespace boost {
namespace numeric {
namespace odeint {


/*
 * ToDo:
 *
 * 1. roll out parameters
 * 2. how to call, with jacobi and system?
 * 3. Introduce adjust size
 * 4. Interfacing for odeint, check if controlled_error_stepper can be used
 * 5. dense output
 */

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

//	void do_step( System &system , Jacobi &jacobi , ... );
//	void do_step( pair< System& , Jacobi& > system , ... );
//	void do_step( pair< System , Jacobi > & , ... );


	template< class System , class Jacobi >
	void do_step( System &system , Jacobi &jacobi , const state_type &x , time_type t , state_type &xout , time_type dt , state_type &xerr )
	{
		const double gamma = 0.25;
		const double d1 = 0.25 , d2 = -0.1043 , d3 = 0.1035 , d4 = 0.3620000000000023e-01;
		const double c2 = 0.386 , c3 = 0.21 , c4 = 0.63;
		const double c21 = -0.5668800000000000e+01;
		const double a21 = 0.1544000000000000e+01;
		const double c31 = -0.2430093356833875e+01 , c32 = -0.2063599157091915e+00;
		const double a31 = 0.9466785280815826e+00 , a32 = 0.2557011698983284e+00;
		const double c41 = -0.1073529058151375e+00 , c42 = -0.9594562251023355e+01 , c43 = -0.2047028614809616e+02;
		const double a41 = 0.3314825187068521e+01 , a42 = 0.2896124015972201e+01 , a43 = 0.9986419139977817e+00;
		const double c51 = 0.7496443313967647e+01 , c52 = -0.1024680431464352e+02 , c53 = -0.3399990352819905e+02 , c54 =  0.1170890893206160e+02;
		const double a51 = 0.1221224509226641e+01 , a52 = 0.6019134481288629e+01 , a53 = 0.1253708332932087e+02 , a54 = -0.6878860361058950e+00 ;
		const double c61 = 0.8083246795921522e+01 , c62 = -0.7981132988064893e+01 , c63 = -0.3152159432874371e+02 , c64 = 0.1631930543123136e+02 , c65 = -0.6058818238834054e+01;

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
    : m_rb4() , m_atol( atol ) , m_rtol( rtol ) ,
      m_first_step( true ) , m_err_old( 0.0 ) , m_dt_old( 0.0 ) ,
      m_last_rejected( false )
	{
	}

	time_type error( const state_type &x , const state_type &xold , const state_type &xerr )
	{
		const size_t n = x.size();
		time_type err = 0.0 , sk = 0.0;
		for( size_t i=0 ; i<n ; ++i )
		{
			sk = m_atol + m_rtol * std::max( std::abs( xold[i] ) , std::abs( x[i] ) );
			err += xerr[i] * xerr[i] / sk / sk;
		}
		return sqrt( err / time_type( n ) );
	}

	time_type last_error( void ) const
	{
		return m_err_old;
	}


	template< class System , class Jacobi >
	boost::numeric::odeint::controlled_step_result
	try_step( System &sys , Jacobi &jacobi , state_type &x , time_type &t , time_type &dt )
	{
		static const time_type safe = 0.9 , fac1 = 5.0 , fac2 = 1.0 / 6.0;

		const size_t n = x.size();
		state_type xnew( n ) , xerr( n );
		m_rb4.do_step( sys , jacobi , x , t , xnew , dt , xerr );
		time_type err = error( xnew , x , xerr );

		time_type fac = std::max( fac2 ,std::min( fac1 , std::pow( err , 0.25 ) / safe ) );
		time_type dt_new = dt / fac;
		if ( err <= 1.0 )
		{
			if( m_first_step )
			{
				m_first_step = false;
			}
			else
			{
				time_type fac_pred = ( m_dt_old / dt ) * pow( err * err / m_err_old , 0.25 ) / safe;
				fac_pred = std::max( fac2 , std::min( fac1 , fac_pred ) );
				fac = std::max( fac , fac_pred );
				dt_new = dt / fac;
			}

			m_dt_old = dt;
			m_err_old = std::max( 0.01 , err );
			if( m_last_rejected )
				dt_new = ( dt >= 0.0 ? std::min( dt_new , dt ) : std::max( dt_new , dt ) );
			t += dt;
			dt = dt_new;
			m_last_rejected = false;
			x = xnew;
			return success_step_size_increased;
		}
		else
		{
			dt = dt_new;
			m_last_rejected = true;
			return step_size_decreased;
		}
	}



private:

	rosenbrock4< time_type > m_rb4;
	time_type m_atol , m_rtol;
	bool m_first_step;
	time_type m_err_old , m_dt_old;
	bool m_last_rejected;
};


}
}
}


#endif /* ROSENBROCK4_HPP_ */
