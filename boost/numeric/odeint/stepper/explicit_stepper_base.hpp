/*
 boost header: BOOST_NUMERIC_ODEINT/explicit_stepper_base.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_EXPLICIT_STEPPER_BASE_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_EXPLICIT_STEPPER_BASE_HPP_INCLUDED

namespace boost {
namespace numeric {
namespace odeint {



struct adjust_size_manually_tag {};
struct adjust_size_initially_tag {};
struct adjust_size_always_tag {};


template<
	class Stepper ,
	unsigned short Order ,
	class State ,
	class Time ,
	class Algebra ,
	class Operations ,
	class AdjustSizePolicy
>
class explicit_stepper_base
{
public:

	typedef State state_type;
	typedef Time time_type;
	typedef Algebra algebra_type;
	typedef Operations operations_type;
	typedef AdjustSizePolicy adjust_size_policy;
	typedef Stepper stepper_type;

	typedef unsigned short order_type;
	static const order_type order_value = Order;

	explicit_stepper_base( void ) : m_is_initialized( false ) { }

    stepper_type& stepper( void ) { return *static_cast< stepper_type* >( this ); }

    const stepper_type& stepper( void ) const {return *static_cast< const stepper_type* >( this );}

    order_type order( void ) const { return order_value; }

	void adjust_size( const state_type &x )
	{
		do_adjust_size( x , typename is_resizeable< state_type >::type() );
	}

protected:

	void do_adjust_size( const state_type &x , boost::true_type )
	{
		this->stepper().adjust_size_impl( x );
	}

	void do_adjust_size( const state_type &x , boost::false_type )
	{
	}


	void adjust_size_by_policy( const state_type &x , adjust_size_manually_tag )
	{
	}

	void adjust_size_by_policy( const state_type &x , adjust_size_initially_tag )
	{
		if( !m_is_initialized )
		{
			do_adjust_size( x , typename is_resizeable< state_type >::type() );
			m_is_initialized = true;
		}
	}

	void adjust_size_by_policy( const state_type &x , adjust_size_always_tag )
	{
		do_adjust_size( x , typename is_resizeable< state_type >::type() );
	}

private:

	bool m_is_initialized;
};



} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_EXPLICIT_STEPPER_BASE_HPP_INCLUDED
