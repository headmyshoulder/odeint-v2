/*
 boost header: numeric/odeint/dense_output_explicit_euler.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_DENSE_OUTPUT_EXPLICIT_EULER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_DENSE_OUTPUT_EXPLICIT_EULER_HPP_INCLUDED

#include <utility>

#include <boost/numeric/odeint/stepper/explicit_euler.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template<
    class State ,
    class Time = double ,
	class Algebra = standard_algebra ,
	class Operations = standard_operations ,
	class AdjustSizePolicy = adjust_size_initially_tag
	>
class dense_output_explicit_euler
{
public:

	typedef State state_type;
	typedef Time time_type;
	typedef Algebra algebra_type;
	typedef Operations operations_type;
	typedef AdjustSizePolicy adjust_size_policy;
	typedef explicit_euler< state_type , time_type , algebra_type , operations_type , adjust_size_policy > stepper_type;

	dense_output_explicit_euler( void )
	: m_euler() , m_size_adjuster() ,
	  m_x1() , m_x2() , m_current_state( &m_x1 ) , m_old_state( &m_x2 ) ,
	  m_t( 0.0 ) , m_t_old( 0.0 ) , m_dt( 1.0 )
	{
		boost::numeric::odeint::construct( m_x1 );
		boost::numeric::odeint::construct( m_x2 );
		m_size_adjuster.register_state( 0 , m_x1 );
		m_size_adjuster.register_state( 1 , m_x2 );
	}

	~dense_output_explicit_euler( void )
	{
		boost::numeric::odeint::destruct( m_x1 );
		boost::numeric::odeint::destruct( m_x2 );
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
		m_euler.do_step( system , *m_current_state , m_t , *m_old_state , m_dt );
		m_t_old = m_t;
		m_t += m_dt;
		std::swap( m_current_state , m_old_state );
		return std::make_pair( m_t_old , m_dt );
	}

	void calc_state( time_type t , state_type &x )
	{
		time_type delta = t - m_t_old;
		algebra_type::for_each3( x , *m_old_state , m_euler.m_dxdt , typename operations_type::template scale_sum2< time_type , time_type >( 1.0 , delta ) );
	}

	void adjust_size( const state_type &x )
	{
		m_size_adjuster.adjust_size( x );
		m_euler.adjust_size( x );
	}


	const state_type& current_state( void ) const { return *m_current_state; }
	const time_type& current_time( void ) const { return m_t; }
	const time_type& previous_state( void ) const { return *m_old_state; }
	const time_type& previous_time( void ) const { return m_t_old; }

private:

	stepper_type m_euler;
	size_adjuster< state_type , 2 > m_size_adjuster;
	state_type m_x1 , m_x2;
	state_type *m_current_state , *m_old_state;
	time_type m_t , m_t_old , m_dt;
};


} // namespace odeint
} // namespace numeric
} // namespace boost


#endif //BOOST_NUMERIC_ODEINT_DENSE_OUTPUT_EXPLICIT_EULER_HPP_INCLUDED
