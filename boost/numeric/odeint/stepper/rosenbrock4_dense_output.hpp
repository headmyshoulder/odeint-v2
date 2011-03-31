/*
 * rosenbrock4_dense_output.hpp
 *
 *  Created on: Mar 29, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_ROSENBROCK4_DENSE_OUTPUT_HPP_
#define STEPPER_ROSENBROCK4_DENSE_OUTPUT_HPP_

#include <utility>

#include <boost/numeric/odeint/stepper/rosenbrock4_controller.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template< class ControlledStepper >
class rosenbrock4_dense_output
{
	void initialize_variables( void )
	{
		m_state_adjuster.register_state( 0 , m_x1 );
		m_state_adjuster.register_state( 1 , m_x2 );
	}

	void copy_variables( const rosenbrock4_dense_output &rb )
	{
		m_stepper = rb.m_stepper;
		m_x1 = rb.m_x1;
		m_x2 = rb.m_x2;
		if( rb.m_current_state == ( & ( rb.m_x1 ) ) )
		{
			m_current_state = &m_x1;
			m_old_state = &m_x2;
		}
		else
		{
			m_current_state = &m_x2;
			m_old_state = &m_x1;
		}
		m_t = rb.m_t;
		m_t_old = rb.m_t_old;
		m_dt = rb.m_dt;
	}

public:

	typedef ControlledStepper controlled_stepper_type;
	typedef typename controlled_stepper_type::stepper_type stepper_type;
	typedef typename stepper_type::value_type value_type;
	typedef typename stepper_type::state_type state_type;
	typedef typename stepper_type::time_type time_type;
	typedef typename stepper_type::deriv_type deriv_type;
	typedef typename stepper_type::adjust_size_policy adjust_size_policy;
	typedef dense_output_stepper_tag stepper_category;

	rosenbrock4_dense_output( const controlled_stepper_type &stepper = controlled_stepper_type() )
	: m_stepper( stepper ) , m_state_adjuster() ,
	  m_x1() , m_x2() , m_current_state( &m_x1 ) , m_old_state( &m_x2 ) ,
	  m_t() , m_t_old() , m_dt()
	{
		initialize_variables();
	}

	rosenbrock4_dense_output( const rosenbrock4_dense_output &rb )
	: m_stepper() , m_state_adjuster() ,
	  m_x1() , m_x2() , m_current_state( &m_x1 ) , m_old_state( &m_x2 ) ,
	  m_t() , m_t_old() , m_dt()
	{
		initialize_variables();
		copy_variables( rb );
	}

	rosenbrock4_dense_output& operator=( const rosenbrock4_dense_output &rb )
	{
		copy_variables( rb );
		return *this;
	}



	template< class StateType >
	void initialize( const StateType &x0 , const time_type &t0 , const time_type &dt0 )
	{
		m_state_adjuster.adjust_size_by_policy( x0 , adjust_size_policy() );
		*m_current_state = x0;
		m_t = t0;
		m_dt = dt0;
	}

	template< class System >
	std::pair< time_type , time_type > do_step( System system )
	{
		const size_t max_count = 1000;

		controlled_step_result res = step_size_decreased;
		m_t_old = m_t;
		size_t count = 0;
		do
		{
			res = m_stepper.try_step( system , *m_current_state , m_t , *m_old_state , m_dt );
			if( count++ == max_count )
				throw std::overflow_error( "dense_output_controlled_explicit_fsal : too much iterations!");
		}
		while( res == step_size_decreased );
		m_stepper.stepper().prepare_dense_output();
		std::swap( m_current_state , m_old_state );
		return std::make_pair( m_t_old , m_t );
	}


	/*
	 * The two overloads are needed in order to solve the forwarding problem.
	 */
	template< class StateOut >
	void calc_state( const time_type &t , StateOut &x )
	{
		m_stepper.stepper().calc_state( t , x , *m_old_state , m_t_old , *m_current_state , m_t );
	}

	template< class StateOut >
	void calc_state( const time_type &t , const StateOut &x )
	{
		m_stepper.stepper().calc_state( t , x , *m_old_state , m_t_old , *m_current_state , m_t );
	}


	template< class StateType >
	void adjust_size( const StateType &x )
	{
		m_stepper.adjust_size( x );
		m_state_adjuster.adjust_size( x );
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
	size_adjuster< state_type , 2 > m_state_adjuster;
	state_type m_x1 , m_x2;
	state_type *m_current_state , *m_old_state;
	time_type m_t , m_t_old , m_dt;
};



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* STEPPER_ROSENBROCK4_DENSE_OUTPUT_HPP_ */
