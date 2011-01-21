/*
 boost header: boost/numeric/odeint/dense_output_dopri5.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_DENSE_OUTPUT_DOPRI5_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_DENSE_OUTPUT_DOPRI5_HPP_INCLUDED

//#include <iostream>
//#define tab "\t"
//using std::cout;
//using std::cerr;
//using std::clog;
//using std::endl;

#include <stdexcept>

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/ref.hpp>

#include <boost/numeric/odeint/stepper/explicit_error_dopri5.hpp>

namespace boost {
namespace numeric {
namespace odeint {


template< class ControlledStepper >
class dense_output_dopri5
{
public:

	typedef ControlledStepper stepper_type;

	typedef typename stepper_type::stepper_type dopri5_type;
	typedef typename dopri5_type::state_type state_type;
	typedef typename dopri5_type::value_type value_type;
	typedef typename dopri5_type::deriv_type deriv_type;
	typedef typename dopri5_type::time_type time_type;
	typedef typename dopri5_type::algebra_type algebra_type;
	typedef typename dopri5_type::operations_type operations_type;
	typedef typename dopri5_type::adjust_size_policy adjust_size_policy;

	BOOST_STATIC_ASSERT(( boost::is_same<
				dopri5_type ,
				explicit_error_dopri5< state_type , value_type , deriv_type , time_type , algebra_type , operations_type , adjust_size_policy >
		>::value ));

	dense_output_dopri5( stepper_type &stepper )
	: m_stepper( stepper ) , m_size_adjuster() ,
	  m_x1() , m_x2() , m_dxdt1() , m_dxdt2() ,
	  m_current_state( &m_x1 ) , m_old_state( &m_x2 ) ,
	  m_current_deriv( &m_dxdt1 ) , m_old_deriv( &m_dxdt2 ) ,
	  m_t( 0.0 ) , m_t_old( 0.0 ) , m_dt( 1.0 ) ,
	  m_is_deriv_initialized( false )
	{
		boost::numeric::odeint::construct( m_x1 );
		boost::numeric::odeint::construct( m_x2 );
		boost::numeric::odeint::construct( m_dxdt1 );
		boost::numeric::odeint::construct( m_dxdt2 );
		m_size_adjuster.register_state( 0 , m_x1 );
		m_size_adjuster.register_state( 1 , m_x2 );
		m_size_adjuster.register_state( 1 , m_dxdt1);
		m_size_adjuster.register_state( 1 , m_dxdt2 );
	}

	~dense_output_dopri5( void )
	{
		boost::numeric::odeint::destruct( m_x1 );
		boost::numeric::odeint::destruct( m_x2 );
		boost::numeric::odeint::destruct( m_dxdt1 );
		boost::numeric::odeint::destruct( m_dxdt2 );
	}

	void adjust_size( const state_type &x )
	{
		m_size_adjuster.adjust_size( x );
		m_stepper.adjust_size( x );
		m_is_deriv_initialized = false;
	}

	void initialize( const state_type &x0 , const time_type t0 , const time_type dt0 )
	{
		boost::numeric::odeint::copy( x0 , *m_current_state );
		m_t = t0;
		m_dt = dt0;
		m_is_deriv_initialized = false;
	}

	template< class System >
	std::pair< time_type , time_type > do_step( System system )
	{
		const size_t max_count = 1000;

		if( !m_is_deriv_initialized )
		{
			typename boost::unwrap_reference< System >::type &sys = system;
			sys( *m_current_state , *m_current_deriv , m_t );
			m_is_deriv_initialized = true;
		}

		controlled_step_result res = step_size_decreased;
		m_t_old = m_t;
		size_t count = 0;
		do
		{
			res = m_stepper.try_step( system , *m_current_state , *m_current_deriv , m_t , *m_old_state , *m_old_deriv , m_dt );
			if( count++ == max_count )
				throw std::overflow_error( "dense_output_dopri5 : too much iterations!");
		}
		while( res == step_size_decreased );
		std::swap( m_current_state , m_old_state );
		std::swap( m_current_deriv , m_old_deriv );
		return std::make_pair( m_t_old , m_t );
	}

	/*
	 * Calculates Dense-Output for Dopri5
	 *
	 * See Hairer, Norsett, Wanner: Solving Ordinary Differential Equations, Nonstiff Problems. I, p.191/192
	 *
	 * y(t+theta) = y(t) + h * sum_i^7 b_i(theta) * k_i
	 *
	 * A = theta^2 * ( 3 - 2 theta )
	 * B = theta^2 * ( theta - 1 )
	 * C = theta^2 * ( theta - 1 )^2
	 * D = theta   * ( theta - 1 )^2
	 *
	 * b_1( theta ) = A * b_1 - C * X1( theta ) + D
	 * b_2( theta ) = 0
	 * b_3( theta ) = A * b_3 + C * X3( theta )
	 * b_4( theta ) = A * b_4 - C * X4( theta )
	 * b_5( theta ) = A * b_5 + C * X5( theta )
	 * b_6( theta ) = A * b_6 - C * X6( theta )
	 * b_7( theta ) = B + C * X7( theta )
	 *
	 * An alternative Method is described in:
	 *
	 * www-m2.ma.tum.de/homepages/simeon/numerik3/kap3.ps
	 */
	void calc_state( time_type t , state_type &x )
	{
		const time_type b1 = static_cast<time_type> ( 35.0 ) / static_cast<time_type>( 384.0 );
		const time_type b3 = static_cast<time_type> ( 500.0 ) / static_cast<time_type>( 1113.0 );
		const time_type b4 = static_cast<time_type> ( 125.0 ) / static_cast<time_type>( 192.0 );
		const time_type b5 = static_cast<time_type> ( -2187.0 ) / static_cast<time_type>( 6784.0 );
		const time_type b6 = static_cast<time_type> ( 11.0 ) / static_cast<time_type>( 84.0 );

		time_type dt = ( m_t - m_t_old );
		time_type theta = ( t - m_t_old ) / dt;
		time_type X1 = static_cast< time_type >( 5.0 ) * ( static_cast< time_type >( 2558722523.0 ) - static_cast< time_type >( 31403016.0 ) * theta ) / static_cast< time_type >( 11282082432.0 );
		time_type X3 = static_cast< time_type >( 100.0 ) * ( static_cast< time_type >( 882725551.0 ) - static_cast< time_type >( 15701508.0 ) * theta ) / static_cast< time_type >( 32700410799.0 );
		time_type X4 = static_cast< time_type >( 25.0 ) * ( static_cast< time_type >( 443332067.0 ) - static_cast< time_type >( 31403016.0 ) * theta ) / static_cast< time_type >( 1880347072.0 ) ;
		time_type X5 = static_cast< time_type >( 32805.0 ) * ( static_cast< time_type >( 23143187.0 ) - static_cast< time_type >( 3489224.0 ) * theta ) / static_cast< time_type >( 199316789632.0 );
		time_type X6 = static_cast< time_type >( 55.0 ) * ( static_cast< time_type >( 29972135.0 ) - static_cast< time_type >( 7076736.0 ) * theta ) / static_cast< time_type >( 822651844.0 );
		time_type X7 = static_cast< time_type >( 10.0 ) * ( static_cast< time_type >( 7414447.0 ) - static_cast< time_type >( 829305.0 ) * theta ) / static_cast< time_type >( 29380423.0 );

		time_type theta_m_1 = theta - static_cast< time_type >( 1.0 );
		time_type theta_sq = theta * theta;
		time_type A = theta_sq * ( static_cast< time_type >( 3.0 ) - static_cast< time_type >( 2.0 ) * theta );
		time_type B = theta_sq * theta_m_1;
		time_type C = theta_sq * theta_m_1 * theta_m_1;
		time_type D = theta * theta_m_1 * theta_m_1;

		time_type b1_theta = A * b1 - C * X1 + D;
		time_type b3_theta = A * b3 + C * X3;
		time_type b4_theta = A * b4 - C * X4;
		time_type b5_theta = A * b5 + C * X5;
		time_type b6_theta = A * b6 - C * X6;
		time_type b7_theta = B + C * X7;

		const state_type &k1 = *m_old_deriv;
		const state_type &k3 = dopri5().m_k3;
		const state_type &k4 = dopri5().m_k4;
		const state_type &k5 = dopri5().m_k5;
		const state_type &k6 = dopri5().m_k6;
		const state_type &k7 = *m_current_deriv;

		typename algebra_type::for_each8()( x , *m_old_state , k1 , k3 , k4 , k5 , k6 , k7 ,
			typename operations_type::template scale_sum7< time_type , time_type , time_type , time_type , time_type , time_type , time_type >( 1.0 , dt * b1_theta , dt * b3_theta , dt * b4_theta , dt * b5_theta , dt * b6_theta , dt * b7_theta ) );
	}

	const state_type& current_state( void ) const
	{
		return *m_current_state;
	}

	const time_type& current_time( void ) const
	{
		return m_t;
	}

	const time_type& previous_state( void ) const
	{
		return *m_old_state;
	}

	const time_type& previous_time( void ) const
	{
		return m_t_old;
	}


private:

	dopri5_type& dopri5( void ) { return m_stepper.stepper(); }
	const dopri5_type& dopri5( void ) const { return m_stepper.stepper(); }



	stepper_type &m_stepper;
	size_adjuster< state_type , 4 > m_size_adjuster;
	state_type m_x1 , m_x2 , m_dxdt1 , m_dxdt2;
	state_type *m_current_state , *m_old_state;
	state_type *m_current_deriv , *m_old_deriv;
	time_type m_t , m_t_old , m_dt;
	bool m_is_deriv_initialized;


};

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif //BOOST_BOOST_NUMERIC_ODEINT_DENSE_OUTPUT_DOPRI5_HPP_INCLUDED
