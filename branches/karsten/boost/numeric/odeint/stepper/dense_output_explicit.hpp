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
#include <boost/numeric/odeint/algebra/default_resize.hpp>

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
		m_stepper = dense_output.m_stepper;
		boost::numeric::odeint::copy( dense_output.m_x1 , m_x1 );
		boost::numeric::odeint::copy( dense_output.m_x2 , m_x2 );
		if( dense_output.m_current_state == (&dense_output.m_x1 ) )
		{
			m_current_state = &m_x1;
			m_old_state = &m_x2;
		}
		else
		{
			m_current_state = &m_x2;
			m_old_state = &m_x1;
		}
		m_t = dense_output.m_t;
		m_t_old = dense_output.m_t_old;
		m_dt = dense_output.m_dt;
	}

public:

	/*
	 * We do not need all typedefs.
	 */
	typedef Stepper stepper_type;
	typedef typename stepper_type::state_type state_type;
	typedef typename stepper_type::value_type value_type;
	typedef typename stepper_type::deriv_type deriv_type;
	typedef typename stepper_type::time_type time_type;
	typedef typename stepper_type::algebra_type algebra_type;
	typedef typename stepper_type::operations_type operations_type;
	typedef typename stepper_type::adjust_size_policy adjust_size_policy;

	dense_output_explicit( const stepper_type &stepper = stepper_type() )
	: m_stepper( stepper ) , m_size_adjuster() ,
	  m_x1() , m_x2() , m_current_state( &m_x1 ) , m_old_state( &m_x2 ) ,
	  m_t( 0.0 ) , m_t_old( 0.0 ) , m_dt( 1.0 )
	{
		initialize_variables();
	}

	~dense_output_explicit( void )
	{
		boost::numeric::odeint::destruct( m_x1 );
		boost::numeric::odeint::destruct( m_x2 );
	}

	dense_output_explicit( const dense_output_explicit &dense_output )
	: m_stepper( dense_output.m_stepper ) , m_size_adjuster() ,
	  m_x1() , m_x2() , m_current_state( &m_x1 ) , m_old_state( &m_x2 ) ,
	  m_t( 0.0 ) , m_t_old( 0.0 ) , m_dt( 1.0 )
	{
		initialize_variables();
		copy_variables( dense_output );
	}

	dense_output_explicit& operator=( const dense_output_explicit &dense_output )
	{
		copy_variables( dense_output );
		return *this;
	}

	template< class StateType >
	void initialize( const StateType &x0 , const time_type &t0 , const time_type &dt0 )
	{
		adjust_size_by_policy( x0 );
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

	template< class StateOut >
	void calc_state( const time_type &t , StateOut &x )
	{
		m_stepper.calc_state( x , t , *m_old_state , m_t_old );
	}

	template< class StateType >
	void adjust_size( const StateType &x )
	{
		m_size_adjuster.adjust_size( x );
		m_stepper.adjust_size( x );
	}

	template< class StateType >
	void adjust_size_by_policy( const StateType &x )
	{
		m_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
		m_stepper.adjust_size_by_policy( x );
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
