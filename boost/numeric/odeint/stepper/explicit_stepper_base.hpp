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


#include <boost/numeric/odeint/algebra/standard_resize.hpp>



namespace boost {
namespace numeric {
namespace odeint {



/*
 * Tags to specify resize behavior of steppers
 */
struct adjust_size_manually_tag {};
struct adjust_size_initially_tag {};
struct adjust_size_always_tag {};




/*
 * Adjust size functionality with policies and resizeability
 */
template< class State , class Stepper >
class size_adjuster
{
public:

	size_adjuster( Stepper &stepper ) : m_is_initialized( false ) , m_stepper( stepper ) { }

	void adjust_size( const State &x )
	{
		adjust_size_by_resizeability( x , typename is_resizeable< State >::type() );
	}

	void adjust_size_by_policy( const State &x , adjust_size_manually_tag )
	{
	}

	void adjust_size_by_policy( const State &x , adjust_size_initially_tag )
	{
		if( !m_is_initialized )
		{
			adjust_size_by_resizeability( x , typename is_resizeable< State >::type() );
			m_is_initialized = true;
		}
	}

	void adjust_size_by_policy( const State &x , adjust_size_always_tag )
	{
		adjust_size_by_resizeability( x , typename is_resizeable< State >::type() );
	}



private:


	void adjust_size_by_resizeability( const State &x , boost::true_type )
	{
		m_stepper.adjust_size_impl( x );
	}

	void adjust_size_by_resizeability( const State &x , boost::false_type )
	{
	}


private :

	bool m_is_initialized;
	Stepper &m_stepper;
};





/*
 * base class for explicit steppers
 */
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

    order_type order( void ) const { return order_value; }


	explicit_stepper_base( void ) : m_size_adjuster( *this ) { }


    stepper_type& stepper( void ) { return *static_cast< stepper_type* >( this ); }
    const stepper_type& stepper( void ) const {return *static_cast< const stepper_type* >( this );}



	template< class System >
	void do_step( System system , state_type &x , time_type t , time_type dt )
	{
		m_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
		system( x , m_dxdt ,t );
		this->stepper().do_step_impl( system , x , m_dxdt , t , dt );
	}


	void adjust_size( const state_type &x )
	{
		m_size_adjuster.adjust_size( x );
	}


private:

	typedef explicit_stepper_base< Stepper , Order , State , Time , Algebra , Operations , AdjustSizePolicy > internal_stepper_base_type;
	typedef size_adjuster< state_type , internal_stepper_base_type > base_size_adjuster_type;
	friend class size_adjuster< state_type , internal_stepper_base_type >;

	void adjust_size_impl( const state_type &x )
	{
		boost::numeric::odeint::adjust_size( x , m_dxdt );
	}


	state_type m_dxdt;
	base_size_adjuster_type m_size_adjuster;
};




} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_EXPLICIT_STEPPER_BASE_HPP_INCLUDED
