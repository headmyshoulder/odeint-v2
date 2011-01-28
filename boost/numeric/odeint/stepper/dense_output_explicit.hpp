/*
 * dense_output_eplicit_stepper.hpp
 *
 *  Created on: Jan 28, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DENSE_OUTPUT_EXPLICIT_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_DENSE_OUTPUT_EXPLICIT_HPP_

#include <utility>

#include <boost/numeric/odeint/stepper/size_adjuster.hpp>

namespace boost {
namespace numeric {
namespace odeint {


template
<
	class Stepper
>
class dense_output_explicit
{
private:

	void initialize_variables( void )
	{
		boost::numeric::odeint::construct( m_x1 );
		boost::numeric::odeint::construct( m_x2 );
		m_size_adjuster.register_state( 0 , m_x1 );
		m_size_adjuster.register_state( 1 , m_x2 );
	}

	void copy_variables( const dense_output_explicit &dense_output )
	{
		/* ToDo : implement */
	}

public:

	/*
	 * ToDo : check which types we really need
	 */
	typedef Stepper stepper_type;
	typedef typename stepper_type::state_type state_type;
	typedef typename stepper_type::value_type value_type;
	typedef typename stepper_type::deriv_type deriv_type;
	typedef typename stepper_type::time_type time_type;
	typedef typename stepper_type::algebra_type algebra_type;
	typedef typename stepper_type::operations_type operations_type;
	typedef typename stepper_type::adjust_size_policy adjust_size_policy;


	dense_output_explicit( const stepper_type &stepper )
	: m_stepper( stepper ) , m_size_adjuster() ,
	  m_x1() , m_x2() , m_current_state( &m_x1 ) , m_old_state( &m_x2 ) ,
	  m_t( 0.0 ) , m_t_old( 0.0 ) , m_dt( 1.0 )
	{
		initialize_variables();
	}

	dense_output_explicit( const dense_output_explicit &dense_ouput )
	{
		initialize_variables();
		copy_variables( dense_output );
	}

	dense_output_explicit& operator=( const dense_output_explicit &dense_output )
	{
		copy_variables( dense_output );
		return *this;
	}

	void initialize( const state_type &x0 , const time_type t0 , const time_type dt0 )
	{
		boost::numeric::odeint::copy( x0 , *m_current_state );
		m_t = t0;
		m_dt = dt0;
	}

	template< class System >
	std::pair< time_type , time_type > do_step( System system )
	{
		m_stepper.do_step( system , *m_current_state , m_t , *m_old_state , m_dt );
		m_t_old = m_t;
		m_t += m_dt;
		std::swap( m_current_state , m_old_state );
		return std::make_pair( m_t_old , m_dt );
	}

	void calc_state( time_type t , state_type &x )
	{
		m_stepper.calc_state( x , *m_old_state , t , m_t_old );
//		time_type delta = t - m_t_old;
//		typename algebra_type::for_each3()( x , *m_old_state , m_euler.m_dxdt , typename operations_type::template scale_sum2< time_type , time_type >( 1.0 , delta ) );
	}

	void adjust_size( const state_type &x )
	{
		m_size_adjuster.adjust_size( x );
		m_stepper.adjust_size( x );
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

	stepper_type m_stepper;
	size_adjuster< state_type , 2 > m_size_adjuster;
	state_type m_x1 , m_x2;
	state_type *m_current_state , *m_old_state;
	time_type m_t , m_t_old , m_dt;

};

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif /* BOOST_NUMERIC_ODEINT_STEPPER_DENSE_OUTPUT_EPLICIT_STEPPER_HPP_ */
