/*
 boost header: numeric/odeint/explicit_stepper_base.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_EXPLICIT_STEPPER_BASE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_EXPLICIT_STEPPER_BASE_HPP_INCLUDED

#include <boost/ref.hpp>

#include <boost/numeric/odeint/stepper/size_adjuster.hpp>
#include <boost/numeric/odeint/algebra/default_resize.hpp>


namespace boost {
namespace numeric {
namespace odeint {

/*
 * base class for explicit steppers
 * models the stepper concept
 */
template<
	class Stepper ,
	unsigned short Order ,
	class State ,
	class Value ,
	class Deriv ,
	class Time ,
	class Algebra ,
	class Operations ,
	class AdjustSizePolicy
>
class explicit_stepper_base
{
public:


	typedef State state_type;
	typedef Value value_type;
	typedef Deriv deriv_type;
	typedef Time time_type;
	typedef Algebra algebra_type;
	typedef Operations operations_type;
	typedef AdjustSizePolicy adjust_size_policy;
	typedef Stepper stepper_type;

	typedef explicit_stepper_base< Stepper , Order , State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > internal_stepper_base_type;

	typedef unsigned short order_type;
	static const order_type order_value = Order;


	order_type order( void ) const
    {
    	return order_value;
    }


	explicit_stepper_base( void ) : m_size_adjuster() , m_dxdt()
	{
		boost::numeric::odeint::construct( m_dxdt );
		m_size_adjuster.register_state( 0 , m_dxdt );
	}

	~explicit_stepper_base( void )
	{
		boost::numeric::odeint::destruct( m_dxdt );
	}

	explicit_stepper_base( const explicit_stepper_base &b ) : m_size_adjuster() , m_dxdt()
	{
		boost::numeric::odeint::construct( m_dxdt );
		m_size_adjuster.register_state( 0 , m_dxdt );
		boost::numeric::odeint::copy( b.m_dxdt , m_dxdt );
	}

	explicit_stepper_base& operator=( const explicit_stepper_base &b )
	{
		boost::numeric::odeint::copy( b.m_dxdt , m_dxdt );
		return *this;
	}





	// do_step( sys , x , t , dt )
	template< class System , class StateInOut >
	void do_step( System system , StateInOut &x , const time_type &t , const time_type &dt )
	{
		typename boost::unwrap_reference< System >::type &sys = system;
		m_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
		sys( x , m_dxdt ,t );
		this->stepper().do_step_impl( system , x , m_dxdt , t , x , dt );
	}

	// do_step( sys , x , dxdt , t , dt )
	template< class System , class StateInOut , class DerivIn >
	void do_step( System system , StateInOut &x , const DerivIn &dxdt , const time_type &t , const time_type &dt )
	{
		this->stepper().do_step_impl( system , x , dxdt , t , x , dt );
	}

	// do_step( sys , in , t , out , dt )
	template< class System , class StateIn , class StateOut >
	void do_step( System system , const StateIn &in , const time_type &t , StateOut &out , const time_type &dt )
	{
		typename boost::unwrap_reference< System >::type &sys = system;
		m_size_adjuster.adjust_size_by_policy( in , adjust_size_policy() );
		sys( in , m_dxdt ,t );
		this->stepper().do_step_impl( system , in , m_dxdt , t , out , dt );
	}

	// do_step( sys , in , dxdt , t , out , dt )
	template< class System , class StateIn , class DerivIn , class StateOut >
	void do_step( System system , const StateIn &in , const DerivIn &dxdt , const time_type &t , StateOut &out , const time_type &dt )
	{
		this->stepper().do_step_impl( system , in , dxdt , t , out , dt );
	}



	template< class StateType >
	void adjust_size( const StateType &x )
	{
		m_size_adjuster.adjust_size( x );
	}

	template< class StateType >
	void adjust_size_by_policy( const StateType &x )
	{
		m_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
	}


protected:

    stepper_type& stepper( void )
    {
    	return *static_cast< stepper_type* >( this );
    }

    const stepper_type& stepper( void ) const
    {
    	return *static_cast< const stepper_type* >( this );
    }


	size_adjuster< deriv_type , 1 > m_size_adjuster;
	deriv_type m_dxdt;
};


} // odeint
} // numeric
} // boost

#endif //BOOST_NUMERIC_ODEINT_EXPLICIT_STEPPER_BASE_HPP_INCLUDED
