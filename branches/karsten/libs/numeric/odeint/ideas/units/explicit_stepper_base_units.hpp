/*
 boost header: numeric/odeint/explicit_stepper_base.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_EXPLICIT_STEPPER_BASE_UNITS_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_EXPLICIT_STEPPER_BASE_UNITS_HPP_INCLUDED

#include <boost/noncopyable.hpp>
#include <boost/ref.hpp>

#include <boost/numeric/odeint/util/size_adjuster.hpp>
#include <boost/numeric/odeint/util/construct.hpp>
#include <boost/numeric/odeint/util/destruct.hpp>
#include <boost/numeric/odeint/util/copy.hpp>


#include "adjust_size_units.hpp"


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
	class Deriv ,
	class Value ,
	class Time ,
	class Algebra ,
	class Operations ,
	class AdjustSizePolicy
>
class explicit_stepper_base_units : boost::noncopyable
{
public:


	typedef Deriv deriv_type;
	typedef Value value_type;
	typedef Time time_type;
	typedef Algebra algebra_type;
	typedef Operations operations_type;
	typedef AdjustSizePolicy adjust_size_policy;
	typedef Stepper stepper_type;

	typedef explicit_stepper_base_units< Stepper , Order , Deriv , Value , Time , Algebra , Operations , AdjustSizePolicy > internal_stepper_base_type;

	typedef unsigned short order_type;
	static const order_type order_value = Order;


	order_type order( void ) const
    {
    	return order_value;
    }


	explicit_stepper_base_units( void ) : m_size_adjuster() , m_dxdt()
	{
		boost::numeric::odeint::construct( m_dxdt );
		m_size_adjuster.register_state( 0 , m_dxdt );
	}

	~explicit_stepper_base_units( void )
	{
		boost::numeric::odeint::destruct( m_dxdt );
	}




	// do_step( sys , x , t , dt )
	template< class System , class State >
	void do_step( System system , State &x , const time_type &t , const time_type &dt )
	{
		typename boost::unwrap_reference< System >::type &sys = system;
		m_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
		sys( x , m_dxdt ,t );
		this->stepper().do_step_impl( system , x , m_dxdt , t , x , dt );
	}

	// do_step( sys , x , dxdt , t , dt )
	template< class System , class State >
	void do_step( System system , State &x , const deriv_type dxdt , const time_type &t , const time_type &dt )
	{
		this->stepper().do_step_impl( system , x , dxdt , t , x , dt );
	}

	// do_step( sys , in , t , out , dt )
	template< class System , class State >
	void do_step( System system , const State &in , const time_type &t , State &out , const time_type &dt )
	{
		typename boost::unwrap_reference< System >::type &sys = system;
		m_size_adjuster.adjust_size_by_policy( in , adjust_size_policy() );
		sys( in , m_dxdt ,t );
		this->stepper().do_step_impl( system , in , m_dxdt , t , out , dt );
	}

	// do_step( sys , in , dxdt , t , out , dt )
	template< class System , class State >
	void do_step( System system , const State &in , const deriv_type &dxdt , const time_type &t , State &out , const time_type &dt )
	{
		this->stepper().do_step_impl( system , in , dxdt , t , out , dt );
	}



	template< class State >
	void adjust_size( const State &x )
	{
		m_size_adjuster.adjust_size( x );
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


	size_adjuster_units< deriv_type , 1 > m_size_adjuster;
	deriv_type m_dxdt;
};


} // odeint
} // numeric
} // boost

#endif //BOOST_NUMERIC_ODEINT_EXPLICIT_STEPPER_BASE_UNITS_HPP_INCLUDED
