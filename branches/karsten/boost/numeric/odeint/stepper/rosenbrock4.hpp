/*
 * rosenbrock4.hpp
 *
 *  Created on: Jan 30, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_ROSENBROCK4_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_ROSENBROCK4_HPP_

#include <boost/ref.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/util/size_adjuster.hpp>

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
 * 6. roll out matrix_vector_adjuster
 */
template< class Value >
class rosenbrock4
{
public:

	typedef Value value_type;
    typedef boost::numeric::ublas::vector< value_type > state_type;
    typedef state_type deriv_type;
    typedef value_type time_type;
    typedef boost::numeric::ublas::matrix< time_type > matrix_type;
    typedef boost::numeric::ublas::permutation_matrix< size_t > pmatrix_type;

	rosenbrock4( void )
	{
	}

	template< class System >
	void do_step( System system , const state_type &x , time_type t , state_type &xout , time_type dt , state_type &xerr )
	{
		const value_type gamma = 0.25;
		const value_type d1 = 0.25 , d2 = -0.1043 , d3 = 0.1035 , d4 = 0.3620000000000023e-01;
		const value_type c2 = 0.386 , c3 = 0.21 , c4 = 0.63;
		const value_type c21 = -0.5668800000000000e+01;
		const value_type a21 = 0.1544000000000000e+01;
		const value_type c31 = -0.2430093356833875e+01 , c32 = -0.2063599157091915e+00;
		const value_type a31 = 0.9466785280815826e+00 , a32 = 0.2557011698983284e+00;
		const value_type c41 = -0.1073529058151375e+00 , c42 = -0.9594562251023355e+01 , c43 = -0.2047028614809616e+02;
		const value_type a41 = 0.3314825187068521e+01 , a42 = 0.2896124015972201e+01 , a43 = 0.9986419139977817e+00;
		const value_type c51 = 0.7496443313967647e+01 , c52 = -0.1024680431464352e+02 , c53 = -0.3399990352819905e+02 , c54 =  0.1170890893206160e+02;
		const value_type a51 = 0.1221224509226641e+01 , a52 = 0.6019134481288629e+01 , a53 = 0.1253708332932087e+02 , a54 = -0.6878860361058950e+00 ;
		const value_type c61 = 0.8083246795921522e+01 , c62 = -0.7981132988064893e+01 , c63 = -0.3152159432874371e+02 , c64 = 0.1631930543123136e+02 , c65 = -0.6058818238834054e+01;

    	typedef typename boost::unwrap_reference< System >::type system_type;
    	typedef typename boost::unwrap_reference< typename system_type::first_type >::type deriv_func_type;
    	typedef typename boost::unwrap_reference< typename system_type::second_type >::type jacobi_func_type;
    	system_type &sys = system;
    	deriv_func_type &deriv_func = sys.first;
    	jacobi_func_type &jacobi_func = sys.second;


		const size_t n = x.size();
		matrix_type jac( n , n );
		pmatrix_type pm( n );
		state_type dfdt( n ) , dxdt( n );
		deriv_func( x , dxdt , t );
		jacobi_func( x , jac , t , dfdt );

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
        deriv_func( xtmp , dxdtnew , t + c2 * dt );
        for( size_t i=0 ; i<n ; ++i )
        	g2[i] = dxdtnew[i] + dt * d2 * dfdt[i] + c21 * g1[i] / dt;
        boost::numeric::ublas::lu_substitute( jac , pm , g2 );


        for( size_t i=0 ; i<n ; ++i )
        	xtmp[i] = x[i] + a31 * g1[i] + a32 * g2[i];
        deriv_func( xtmp , dxdtnew , t + c3 * dt );
        for( size_t i=0 ; i<n ; ++i )
        	g3[i] = dxdtnew[i] + dt * d3 * dfdt[i] + ( c31 * g1[i] + c32 * g2[i] ) / dt;
        boost::numeric::ublas::lu_substitute( jac , pm , g3 );


        for( size_t i=0 ; i<n ; ++i )
        	xtmp[i] = x[i] + a41 * g1[i] + a42 * g2[i] + a43 * g3[i];
        deriv_func( xtmp , dxdtnew , t + c4 * dt );
        for( size_t i=0 ; i<n ; ++i )
        	g4[i] = dxdtnew[i] + dt * d4 * dfdt[i] + ( c41 * g1[i] + c42 * g2[i] + c43 * g3[i] ) / dt;
        boost::numeric::ublas::lu_substitute( jac , pm , g4 );


        for( size_t i=0 ; i<n ; ++i )
        	xtmp[i] = x[i] + a51 * g1[i] + a52 * g2[i] + a53 * g3[i] + a54 * g4[i];
        deriv_func( xtmp , dxdtnew , t + dt );
        for( size_t i=0 ; i<n ; ++i )
        	g5[i] = dxdtnew[i] + ( c51 * g1[i] + c52 * g2[i] + c53 * g3[i] + c54 * g4[i] ) / dt;
        boost::numeric::ublas::lu_substitute( jac , pm , g5 );

        for( size_t i=0 ; i<n ; ++i )
        	xtmp[i] += g5[i];
        deriv_func( xtmp , dxdtnew , t + dt );
        for( size_t i=0 ; i<n ; ++i )
        	xerr[i] = dxdtnew[i] + ( c61 * g1[i] + c62 * g2[i] + c63 * g3[i] + c64 * g4[i] + c65 * g5[i] ) / dt;
        boost::numeric::ublas::lu_substitute( jac , pm , xerr );

        for( size_t i=0 ; i<n ; ++i )
        	xout[i] = xtmp[i] + xerr[i];
	}

	template< class System >
	void do_step( System system , state_type &x , time_type t , time_type dt , state_type &xerr )
	{
		do_step( system , x , t , x , dt , xerr );
	}


private:

	size_adjuster< state_type , 9 > m_state_adjuster;
    size_adjuster< matrix_type , 1 , matrix_vector_adjuster> m_matrix_adjuster;
    size_adjuster< pmatrix_type , 1 > m_pmatrix_adjuster;

	matrix_type m_jac;
	pmatrix_type m_pm;
	state_type m_dfdt , m_dxdt;
	state_type m_g1 , m_g2 , m_g3 , m_g4 , m_g5;
	state_type m_xtmp , m_dxdtnew;


};


} // namespace odeint
} // namespace numeric
} // namespace boost

#endif /* BOOST_NUMERIC_ODEINT_STEPPER_ROSENBROCK4_HPP_ */
