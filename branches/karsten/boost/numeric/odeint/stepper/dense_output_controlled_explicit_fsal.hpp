/*
 * dense_output_controlled_explicit_fsal.hpp
 *
 *  Created on: Jan 28, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DENSE_OUTPUT_CONTROLLED_EXPLICIT_FSAL_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_DENSE_OUTPUT_CONTROLLED_EXPLICIT_FSAL_HPP_

#include <utility>
#include <stdexcept>

#include <boost/numeric/odeint/stepper/size_adjuster.hpp>

namespace boost {
namespace numeric {
namespace odeint {


template
<
	class ControlledStepper
>
class dense_output_controlled_explicit_fsal
{
private:

	void initialize_variables( void )
	{
		boost::numeric::odeint::construct( m_x1 );
		boost::numeric::odeint::construct( m_x2 );
		boost::numeric::odeint::construct( m_dxdt1 );
		boost::numeric::odeint::construct( m_dxdt2 );
		m_state_adjuster.register_state( 0 , m_x1 );
		m_state_adjuster.register_state( 1 , m_x2 );
		m_deriv_adjuster.register_state( 0 , m_x1 );
		m_deriv_adjuster.register_state( 1 , m_x2 );
	}

	void copy_variables( const dense_output_controlled_explicit_fsal &dense_output )
	{
		m_stepper = dense_output.m_stepper;
		boost::numeric::odeint::copy( dense_output.m_x1 , m_x1 );
		boost::numeric::odeint::copy( dense_output.m_x2 , m_x2 );
		boost::numeric::odeint::copy( dense_output.m_dxdt1 , m_dxdt1 );
		boost::numeric::odeint::copy( dense_output.m_dxdt2 , m_dxdt2 );
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
		if( dense_output.m_current_deriv == ( &dense_output.m_dxdt1 ) )
		{
			m_current_deriv = &m_dxdt1;
			m_old_deriv = &m_dxdt2;
		}
		else
		{
			m_current_deriv = &m_dxdt2;
			m_old_deriv = &m_dxdt1;
		}
		m_t = dense_output.m_t;
		m_t_old = dense_output.m_t_old;
		m_dt = dense_output.m_dt;
		m_is_deriv_initialized = dense_output.m_is_deriv_initialized;
	}

public:

	/*
	 * We do not need all typedefs.
	 */
	typedef ControlledStepper controlled_stepper_type;

	typedef typename controlled_stepper_type::stepper_type stepper_type;
	typedef typename stepper_type::state_type state_type;
	typedef typename stepper_type::value_type value_type;
	typedef typename stepper_type::deriv_type deriv_type;
	typedef typename stepper_type::time_type time_type;
	typedef typename stepper_type::algebra_type algebra_type;
	typedef typename stepper_type::operations_type operations_type;
	typedef typename stepper_type::adjust_size_policy adjust_size_policy;


	dense_output_controlled_explicit_fsal( const controlled_stepper_type &stepper = controlled_stepper_type() )
	: m_stepper( stepper ) ,
	  m_state_adjuster() , m_deriv_adjuster() ,
	  m_x1() , m_x2() , m_current_state( &m_x1 ) , m_old_state( &m_x2 ) ,
	  m_dxdt1() , m_dxdt2() , m_current_deriv( &m_dxdt1 ) , m_old_deriv( &m_dxdt2 ) ,
	  m_t( 0.0 ) , m_t_old( 0.0 ) , m_dt( 1.0 ) , m_is_deriv_initialized( false )
	{
		initialize_variables();
	}

	dense_output_controlled_explicit_fsal( const dense_output_controlled_explicit_fsal &dense_output )
	: m_stepper( dense_output.m_stepper ) ,
	  m_state_adjuster() , m_deriv_adjuster() ,
	  m_x1() , m_x2() , m_current_state( &m_x1 ) , m_old_state( &m_x2 ) ,
	  m_dxdt1() , m_dxdt2() , m_current_deriv( &m_dxdt1 ) , m_old_deriv( &m_dxdt2 ) ,
	  m_t( 0.0 ) , m_t_old( 0.0 ) , m_dt( 1.0 ) , m_is_deriv_initialized( false )
	{
		initialize_variables();
		copy_variables( dense_output );
	}

	~dense_output_controlled_explicit_fsal( void )
	{
		boost::numeric::odeint::destruct( m_x1 );
		boost::numeric::odeint::destruct( m_x2 );
		boost::numeric::odeint::destruct( m_dxdt1 );
		boost::numeric::odeint::destruct( m_dxdt2 );
	}

	dense_output_controlled_explicit_fsal& operator=( const dense_output_controlled_explicit_fsal &dense_output )
	{
		copy_variables( dense_output );
		return *this;
	}

	template< class StateType >
	void initialize( const StateType &x0 , const time_type &t0 , const time_type &dt0 )
	{
		m_state_adjuster.adjust_size_by_policy( x0 , adjust_size_policy() );
		m_deriv_adjuster.adjust_size_by_policy( x0 , adjust_size_policy() );
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
				throw std::overflow_error( "dense_output_controlled_explicit_fsal : too much iterations!");
		}
		while( res == step_size_decreased );
		std::swap( m_current_state , m_old_state );
		std::swap( m_current_deriv , m_old_deriv );
		return std::make_pair( m_t_old , m_t );
	}

	template< class StateOut >
	void calc_state( const time_type &t , StateOut &x )
	{
		m_stepper.stepper().calc_state( t , x , *m_old_state , *m_old_deriv , m_t_old , *m_current_state , *m_current_deriv , m_t );
	}

	template< class StateType >
	void adjust_size( const StateType &x )
	{
		m_state_adjuster.adjust_size( x );
		m_deriv_adjuster.adjust_size( x );
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

	controlled_stepper_type m_stepper;
	size_adjuster< state_type , 2 > m_state_adjuster;
	size_adjuster< deriv_type , 2 > m_deriv_adjuster;
	state_type m_x1 , m_x2;
	state_type *m_current_state , *m_old_state;
	deriv_type m_dxdt1 , m_dxdt2;
	state_type *m_current_deriv , *m_old_deriv;
	time_type m_t , m_t_old , m_dt;
	bool m_is_deriv_initialized;

};

} // namespace odeint
} // namespace numeric
} // namespace boost



#endif /* BOOST_NUMERIC_ODEINT_STEPPER_DENSE_OUTPUT_CONTROLLED_EXPLICIT_FSAL_HPP_ */
