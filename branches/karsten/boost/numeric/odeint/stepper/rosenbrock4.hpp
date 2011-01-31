/*
 * rosenbrock4.hpp
 *
 *  Created on: Jan 30, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_ROSENBROCK4_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_ROSENBROCK4_HPP_

#include <boost/ref.hpp>
#include <boost/noncopyable.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/util/size_adjuster.hpp>
#include <boost/numeric/odeint/util/matrix_vector_adjust_size.hpp>
#include <boost/numeric/odeint/util/ublas_resize.hpp>
#include <boost/numeric/odeint/util/ublas_permutation_matrix_resize.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>


namespace boost {
namespace numeric {
namespace odeint {


/*
 * ToDo:
 *
 * 2. Interfacing for odeint, check if controlled_error_stepper can be used
 * 3. dense output
 */



template< class Value >
struct default_rosenbrock_coefficients : boost::noncopyable
{
	typedef Value value_type;

	default_rosenbrock_coefficients( void )
	: gamma ( 0.25 ) ,
	  d1 ( 0.25 ) , d2 ( -0.1043 ) , d3 ( 0.1035 ) , d4 ( 0.3620000000000023e-01 ) ,
	  c2 ( 0.386 ) , c3 ( 0.21 ) , c4 ( 0.63 ) ,
	  c21 ( -0.5668800000000000e+01 ) ,
	  a21 ( 0.1544000000000000e+01 ) ,
	  c31 ( -0.2430093356833875e+01 ) , c32 ( -0.2063599157091915e+00 ) ,
	  a31 ( 0.9466785280815826e+00 ) , a32 ( 0.2557011698983284e+00 ) ,
	  c41 ( -0.1073529058151375e+00 ) , c42 ( -0.9594562251023355e+01 ) , c43 ( -0.2047028614809616e+02 ) ,
	  a41 ( 0.3314825187068521e+01 ) , a42 ( 0.2896124015972201e+01 ) , a43 ( 0.9986419139977817e+00 ) ,
	  c51 ( 0.7496443313967647e+01 ) , c52 ( -0.1024680431464352e+02 ) , c53 ( -0.3399990352819905e+02 ) , c54 (  0.1170890893206160e+02 ) ,
	  a51 ( 0.1221224509226641e+01 ) , a52 ( 0.6019134481288629e+01 ) , a53 ( 0.1253708332932087e+02 ) , a54 ( -0.6878860361058950e+00 ) ,
	  c61 ( 0.8083246795921522e+01 ) , c62 ( -0.7981132988064893e+01 ) , c63 ( -0.3152159432874371e+02 ) , c64 ( 0.1631930543123136e+02 ) , c65 ( -0.6058818238834054e+01 )
	  {}

	const value_type gamma;
	const value_type d1 , d2 , d3 , d4;
	const value_type c2 , c3 , c4;
	const value_type c21 ;
	const value_type a21;
	const value_type c31 , c32;
	const value_type a31 , a32;
	const value_type c41 , c42 , c43;
	const value_type a41 , a42 , a43;
	const value_type c51 , c52 , c53 , c54;
	const value_type a51 , a52 , a53 , a54;
	const value_type c61 , c62 , c63 , c64 , c65;
};



template< class Value , class Coefficients = default_rosenbrock_coefficients< Value > , class AdjustSizePolicy = adjust_size_initially_tag >
class rosenbrock4
{
private:

	void initialize( void )
	{
		m_matrix_adjuster.register_state( 0 , m_jac );

		m_pmatrix_adjuster.register_state( 0 , m_pm );

		m_state_adjuster.register_state( 0 , m_dfdt );
		m_state_adjuster.register_state( 1 , m_dxdt );
		m_state_adjuster.register_state( 2 , m_g1 );
		m_state_adjuster.register_state( 3 , m_g2 );
		m_state_adjuster.register_state( 4 , m_g3 );
		m_state_adjuster.register_state( 5 , m_g4 );
		m_state_adjuster.register_state( 6 , m_g5 );
		m_state_adjuster.register_state( 7 , m_xtmp );
		m_state_adjuster.register_state( 8 , m_dxdtnew );
	}

	void copy( const rosenbrock4 &rb )
	{
		m_jac = rb.m_jac;
		m_pm = rb.m_pm;
		m_dfdt =rb.m_dfdt;
		m_dxdt = rb.m_dxdt;
		m_g1 = rb.m_g1;
		m_g2 = rb.m_g2;
		m_g3 = rb.m_g3;
		m_g4 = rb.m_g4;
		m_g5 = rb.m_g5;
		m_xtmp = rb.m_xtmp;
		m_dxdtnew = rb.m_dxdtnew;
	}

public:

	typedef Value value_type;
    typedef boost::numeric::ublas::vector< value_type > state_type;
    typedef state_type deriv_type;
    typedef value_type time_type;
    typedef boost::numeric::ublas::matrix< value_type > matrix_type;
    typedef boost::numeric::ublas::permutation_matrix< size_t > pmatrix_type;
    typedef AdjustSizePolicy adjust_size_policy;
    typedef Coefficients rosenbrock_coefficients;


	rosenbrock4( void )
	: m_state_adjuster() , m_matrix_adjuster() , m_pmatrix_adjuster() ,
	  m_jac() , m_pm( 1 ) ,
	  m_dfdt() , m_dxdt() ,
	  m_g1() , m_g2() , m_g3() , m_g4() , m_g5() ,
	  m_xtmp() , m_dxdtnew() ,
	  m_coef()
    {
		initialize();
	}

	rosenbrock4( const rosenbrock4 &rb )
	: m_state_adjuster() , m_matrix_adjuster() , m_pmatrix_adjuster() ,
	  m_jac() , m_pm( 1 ) ,
	  m_dfdt() , m_dxdt() ,
	  m_g1() , m_g2() , m_g3() , m_g4() , m_g5() ,
	  m_xtmp() , m_dxdtnew() ,
	  m_coef()
	{
		initialize();
		copy( rb );
	}

	rosenbrock4& operator=( const rosenbrock4 &rb )
	{
		copy( rb );
	}

	template< class System >
	void do_step( System system , const state_type &x , time_type t , state_type &xout , time_type dt , state_type &xerr )
	{
		// get the systen and jacobi function
    	typedef typename boost::unwrap_reference< System >::type system_type;
    	typedef typename boost::unwrap_reference< typename system_type::first_type >::type deriv_func_type;
    	typedef typename boost::unwrap_reference< typename system_type::second_type >::type jacobi_func_type;
    	system_type &sys = system;
    	deriv_func_type &deriv_func = sys.first;
    	jacobi_func_type &jacobi_func = sys.second;

    	const size_t n = x.size();

    	// adjust size
		m_matrix_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
		m_pmatrix_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
		for( size_t i=0 ; i<n ; ++i )
			m_pm( i ) = i;
		m_state_adjuster.adjust_size_by_policy( x , adjust_size_policy() );


		deriv_func( x , m_dxdt , t );
		jacobi_func( x , m_jac , t , m_dfdt );

		m_jac *= -1.0;
		m_jac += 1.0 / m_coef.gamma / dt * boost::numeric::ublas::identity_matrix< value_type >( n );
        boost::numeric::ublas::lu_factorize( m_jac , m_pm );

        for( size_t i=0 ; i<n ; ++i )
        	m_g1[i] = m_dxdt[i] + dt * m_coef.d1 * m_dfdt[i];
        boost::numeric::ublas::lu_substitute( m_jac , m_pm , m_g1 );


        for( size_t i=0 ; i<n ; ++i )
        	m_xtmp[i] = x[i] + m_coef.a21 * m_g1[i];
        deriv_func( m_xtmp , m_dxdtnew , t + m_coef.c2 * dt );
        for( size_t i=0 ; i<n ; ++i )
        	m_g2[i] = m_dxdtnew[i] + dt * m_coef.d2 * m_dfdt[i] + m_coef.c21 * m_g1[i] / dt;
        boost::numeric::ublas::lu_substitute( m_jac , m_pm , m_g2 );


        for( size_t i=0 ; i<n ; ++i )
        	m_xtmp[i] = x[i] + m_coef.a31 * m_g1[i] + m_coef.a32 * m_g2[i];
        deriv_func( m_xtmp , m_dxdtnew , t + m_coef.c3 * dt );
        for( size_t i=0 ; i<n ; ++i )
        	m_g3[i] = m_dxdtnew[i] + dt * m_coef.d3 * m_dfdt[i] + ( m_coef.c31 * m_g1[i] + m_coef.c32 * m_g2[i] ) / dt;
        boost::numeric::ublas::lu_substitute( m_jac , m_pm , m_g3 );


        for( size_t i=0 ; i<n ; ++i )
        	m_xtmp[i] = x[i] + m_coef.a41 * m_g1[i] + m_coef.a42 * m_g2[i] + m_coef.a43 * m_g3[i];
        deriv_func( m_xtmp , m_dxdtnew , t + m_coef.c4 * dt );
        for( size_t i=0 ; i<n ; ++i )
        	m_g4[i] = m_dxdtnew[i] + dt * m_coef.d4 * m_dfdt[i] + ( m_coef.c41 * m_g1[i] + m_coef.c42 * m_g2[i] + m_coef.c43 * m_g3[i] ) / dt;
        boost::numeric::ublas::lu_substitute( m_jac , m_pm , m_g4 );


        for( size_t i=0 ; i<n ; ++i )
        	m_xtmp[i] = x[i] + m_coef.a51 * m_g1[i] + m_coef.a52 * m_g2[i] + m_coef.a53 * m_g3[i] + m_coef.a54 * m_g4[i];
        deriv_func( m_xtmp , m_dxdtnew , t + dt );
        for( size_t i=0 ; i<n ; ++i )
        	m_g5[i] = m_dxdtnew[i] + ( m_coef.c51 * m_g1[i] + m_coef.c52 * m_g2[i] + m_coef.c53 * m_g3[i] + m_coef.c54 * m_g4[i] ) / dt;
        boost::numeric::ublas::lu_substitute( m_jac , m_pm , m_g5 );

        for( size_t i=0 ; i<n ; ++i )
        	m_xtmp[i] += m_g5[i];
        deriv_func( m_xtmp , m_dxdtnew , t + dt );
        for( size_t i=0 ; i<n ; ++i )
        	xerr[i] = m_dxdtnew[i] + ( m_coef.c61 * m_g1[i] + m_coef.c62 * m_g2[i] + m_coef.c63 * m_g3[i] + m_coef.c64 * m_g4[i] + m_coef.c65 * m_g5[i] ) / dt;
        boost::numeric::ublas::lu_substitute( m_jac , m_pm , xerr );

        for( size_t i=0 ; i<n ; ++i )
        	xout[i] = m_xtmp[i] + xerr[i];
	}

	template< class System >
	void do_step( System system , state_type &x , time_type t , time_type dt , state_type &xerr )
	{
		do_step( system , x , t , x , dt , xerr );
	}



	template< class StateType >
	void adjust_size( const StateType &x )
	{
		m_state_adjuster.adjust_size( x );
		m_matrix_adjuster.adjust_size( x );
		m_pmatrix_adjuster.adjust_size( x );
	}


private:

	size_adjuster< state_type , 9 > m_state_adjuster;
    size_adjuster< matrix_type , 1 , matrix_vector_adjust_size > m_matrix_adjuster;
    size_adjuster< pmatrix_type , 1 > m_pmatrix_adjuster;

	matrix_type m_jac;
	pmatrix_type m_pm;
	state_type m_dfdt , m_dxdt;
	state_type m_g1 , m_g2 , m_g3 , m_g4 , m_g5;
	state_type m_xtmp , m_dxdtnew;

	rosenbrock_coefficients m_coef;
};


} // namespace odeint
} // namespace numeric
} // namespace boost

#endif /* BOOST_NUMERIC_ODEINT_STEPPER_ROSENBROCK4_HPP_ */
