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

#include <boost/ref.hpp>
#include <boost/bind.hpp>

#include <boost/numeric/odeint/util/copy.hpp>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

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

    void copy_pointers( const dense_output_controlled_explicit_fsal &dense_output )
	{
		if( dense_output.m_current_state == (&dense_output.m_x1.m_v ) )
		{
			m_current_state = &m_x1.m_v;
			m_old_state = &m_x2.m_v;
		}
		else
		{
			m_current_state = &m_x2.m_v;
			m_old_state = &m_x1.m_v;
		}
		if( dense_output.m_current_deriv == ( &dense_output.m_dxdt1.m_v ) )
		{
			m_current_deriv = &m_dxdt1.m_v;
			m_old_deriv = &m_dxdt2.m_v;
		}
		else
		{
			m_current_deriv = &m_dxdt2.m_v;
			m_old_deriv = &m_dxdt1.m_v;
		}
	}

public:

	/*
	 * We do not need all typedefs.
	 */
	typedef ControlledStepper controlled_stepper_type;

	typedef typename controlled_stepper_type::stepper_type stepper_type;
	typedef typename stepper_type::state_type state_type;
	typedef typename stepper_type::wrapped_state_type wrapped_state_type;
	typedef typename stepper_type::value_type value_type;
	typedef typename stepper_type::deriv_type deriv_type;
	typedef typename stepper_type::wrapped_deriv_type wrapped_deriv_type;
	typedef typename stepper_type::time_type time_type;
	typedef typename stepper_type::algebra_type algebra_type;
	typedef typename stepper_type::operations_type operations_type;
	typedef typename stepper_type::resizer_type resizer_type;
	typedef dense_output_stepper_tag stepper_category;
	typedef dense_output_controlled_explicit_fsal< ControlledStepper > dense_output_stepper_type;

	dense_output_controlled_explicit_fsal( const controlled_stepper_type &stepper = controlled_stepper_type() )
	: m_stepper( stepper ) ,
	  m_current_state( &m_x1.m_v ) , m_old_state( &m_x2.m_v ) ,
	  m_current_deriv( &m_dxdt1.m_v ) , m_old_deriv( &m_dxdt2.m_v ) ,
	  m_is_deriv_initialized( false )
	{ }

	dense_output_controlled_explicit_fsal( const dense_output_controlled_explicit_fsal &dense_output )
	: m_stepper( dense_output.m_stepper ) , 
      m_x1( dense_output.m_x1 ) , m_x2( dense_output.m_x2 ) ,
      m_dxdt1( dense_output.m_dxdt1 ) , m_dxdt2( dense_output.m_dxdt2 ) ,
	  m_current_state( &m_x1.m_v ) , m_old_state( &m_x2.m_v ) ,
	  m_current_deriv( &m_dxdt1.m_v ) , m_old_deriv( &m_dxdt2.m_v ) ,
      m_t( dense_output.m_t ) , m_t_old( dense_output.m_t_old ) , m_dt( dense_output.m_dt ) , m_is_deriv_initialized( dense_output.m_is_deriv_initialized )
	{
		copy_pointers( dense_output );
	}

	dense_output_controlled_explicit_fsal& operator=( const dense_output_controlled_explicit_fsal &dense_output )
	{
		copy_pointers( dense_output );
		return *this;
	}

	template< class StateType >
	void initialize( const StateType &x0 , const time_type &t0 , const time_type &dt0 )
	{
	    m_resizer.adjust_size( x0 , boost::bind( &dense_output_stepper_type::resize< StateType > , boost::ref( *this ) , _1 ) );
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


	/*
	 * The two overloads are needed in order to solve the forwarding problem.
	 */
	template< class StateOut >
	void calc_state( const time_type &t , StateOut &x )
	{
		m_stepper.stepper().calc_state( t , x , *m_old_state , *m_old_deriv , m_t_old , *m_current_state , *m_current_deriv , m_t );
	}

	template< class StateOut >
	void calc_state( const time_type &t , const StateOut &x )
	{
		m_stepper.stepper().calc_state( t , x , *m_old_state , *m_old_deriv , m_t_old , *m_current_state , *m_current_deriv , m_t );
	}


	template< class StateIn >
	bool resize( const StateIn &x )
	{
	    bool resized = false;
	    resized |= adjust_size_by_resizeability( m_x1 , x , typename wrapped_state_type::is_resizeable() );
	    resized |= adjust_size_by_resizeability( m_x2 , x , typename wrapped_state_type::is_resizeable() );
	    resized |= adjust_size_by_resizeability( m_dxdt1 , x , typename wrapped_deriv_type::is_resizeable() );
	    resized |= adjust_size_by_resizeability( m_dxdt2 , x , typename wrapped_deriv_type::is_resizeable() );
	    return resized;
	}


	template< class StateType >
	void adjust_size( const StateType &x )
	{
	    resize( x );
	    m_stepper.stepper().resize( x );
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

	const time_type& current_time_step( void ) const
	{
		return m_dt;
	}


private:

	controlled_stepper_type m_stepper;
	resizer_type m_resizer;
	wrapped_state_type m_x1 , m_x2;
	state_type *m_current_state , *m_old_state;
	wrapped_deriv_type m_dxdt1 , m_dxdt2;
	deriv_type *m_current_deriv , *m_old_deriv;
	time_type m_t , m_t_old , m_dt;
	bool m_is_deriv_initialized;

};

} // namespace odeint
} // namespace numeric
} // namespace boost



#endif /* BOOST_NUMERIC_ODEINT_STEPPER_DENSE_OUTPUT_CONTROLLED_EXPLICIT_FSAL_HPP_ */
