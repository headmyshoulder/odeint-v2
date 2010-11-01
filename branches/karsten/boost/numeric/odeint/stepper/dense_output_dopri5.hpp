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

#include <iostream>
#define tab "\t"
using std::cout;
using std::cerr;
using std::clog;
using std::endl;

#include <stdexcept>

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

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
	typedef typename dopri5_type::time_type time_type;
	typedef typename dopri5_type::algebra_type algebra_type;
	typedef typename dopri5_type::operations_type operations_type;
	typedef typename dopri5_type::adjust_size_policy adjust_size_policy;

	BOOST_STATIC_ASSERT(( boost::is_same<
				dopri5_type ,
				explicit_error_dopri5< state_type , time_type , algebra_type , operations_type , adjust_size_policy >
		>::value ));

	dense_output_dopri5( stepper_type &stepper )
	: m_stepper( stepper ) , m_size_adjuster() ,
	  m_x1() , m_x2() , m_current_state( &m_x1 ) , m_old_state( &m_x2 ) ,
	  m_t( 0.0 ) , m_t_old( 0.0 ) , m_dt( 1.0 )
	{
		boost::numeric::odeint::construct( m_x1 );
		boost::numeric::odeint::construct( m_x2 );
		m_size_adjuster.register_state( 0 , m_x1 );
		m_size_adjuster.register_state( 1 , m_x2 );
	}

	~dense_output_dopri5( void )
	{
		boost::numeric::odeint::destruct( m_x1 );
		boost::numeric::odeint::destruct( m_x2 );
	}

	void adjust_size( const state_type &x )
	{
		m_size_adjuster.adjust_size( x );
		m_stepper.adjust_size( x );
	}

	void initialize( const state_type &x0 , const time_type t0 , const time_type dt0 )
	{
		boost::numeric::odeint::copy( x0 , *m_current_state );
		m_t = t0;
		m_dt = dt0;
	}

	template< class System >
	std::pair< time_type , time_type > do_step( System &system )
	{
		const size_t max_count = 1000;
		controlled_step_result res;
		m_t_old = m_t;
		size_t count = 0;
		do
		{
			res = m_stepper.try_step( system , *m_current_state , m_t , *m_old_state , m_dt );
			if( count++ == max_count )
				throw std::overflow_error( "dense_output_dopri5 : too much iterations!");
		}
		while( res == step_size_decreased );
		std::swap( m_current_state , m_old_state );
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
	 * b_1( theta ) = A * b_1 + C * X1( theta ) + D
	 * b_2( theta ) = 0
	 * b_3( theta ) = A * b_3 + C * X3( theta )
	 * b_4( theta ) = A * b_4 + C * X4( theta )
	 * b_5( theta ) = A * b_5 + C * X5( theta )
	 * b_6( theta ) = A * b_6 + C * X6( theta )
	 * b_7( theta ) = B + C * X7( theta )
	 */
	void calc_state( time_type t , state_type &x )
	{

		const time_type b1 = static_cast<time_type> ( 35.0 ) / static_cast<time_type>( 384.0 );
		const time_type b3 = static_cast<time_type> ( 500.0 ) / static_cast<time_type>( 1113.0 );
		const time_type b4 = static_cast<time_type> ( 125.0 ) / static_cast<time_type>( 192.0 );
		const time_type b5 = static_cast<time_type> ( -2187.0 ) / static_cast<time_type>( 6784.0 );
		const time_type b6 = static_cast<time_type> ( 11.0 ) / static_cast<time_type>( 84.0 );


		time_type theta = t - m_t_old;
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

//		time_type b1_theta = A *





//		algebra_type::for_each3( x , *m_old_state , m_euler.m_dxdt , typename operations_type::scale_sum2( 1.0 , delta ) );
	}

	const state_type& current_state( void ) const { return *m_current_state; }
	const time_type& current_time( void ) const { return m_t; }
	const time_type& previous_state( void ) const { return *m_old_state; }
	const time_type& previous_time( void ) const { return m_t_old; }


private:



	stepper_type &m_stepper;
	size_adjuster< state_type , 2 > m_size_adjuster;
	state_type m_x1 , m_x2;
	state_type *m_current_state , *m_old_state;
	time_type m_t , m_t_old , m_dt;


};

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif //BOOST_BOOST_NUMERIC_ODEINT_DENSE_OUTPUT_DOPRI5_HPP_INCLUDED
