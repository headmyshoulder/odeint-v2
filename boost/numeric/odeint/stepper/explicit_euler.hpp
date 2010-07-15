/*
 boost header: BOOST_NUMERIC_ODEINT/explicit_euler.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_EXPLICIT_EULER_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_EXPLICIT_EULER_HPP_INCLUDED

#include <boost/numeric/odeint/algebra/standard_algebra.hpp>
#include <boost/numeric/odeint/algebra/standard_operations.hpp>

#include <boost/numeric/odeint/stepper/explicit_stepper_base.hpp>
#include <boost/numeric/odeint/stepper/detail/macros.hpp>

namespace boost {
namespace numeric {
namespace odeint {


template<
    class State ,
    class Time = double ,
	class Algebra = standard_algebra< State > ,
	class Operations = standard_operations< Time > ,
	class AdjustSizePolicy = adjust_size_initially_tag
	>
class explicit_euler
: public explicit_stepper_base<
	  explicit_euler< State , Time , Algebra , Operations , AdjustSizePolicy > ,
	  1 , State , Time , Algebra , Operations , AdjustSizePolicy >
{
public :


	BOOST_ODEINT_EXPLICIT_STEPPERS_TYPEDEFS( explicit_euler , 1 );

	template< class System >
	void do_step_impl( System system , state_type &x , const state_type &dxdt , time_type t , time_type dt )
	{
		algebra_type::for_each2( x , dxdt , typename operations_type::increment( dt ) );
	}


//	explicit_euler( void ) : m_size_adjuster( *this ) { }
//
//	void adjust_size( const state_type &x )
//	{
//		m_size_adjuster.adjust_size( x );
//		stepper_base_type::adjust_size( x );
//	}
//
//
//private:
//
//	typedef explicit_euler< State , Time , Algebra , Operations , AdjustSizePolicy > stepper_type;
//	typedef size_adjuster< state_type , stepper_type > size_adjuster_type;
//	friend class size_adjuster< state_type , stepper_type >;
//
//	void adjust_size_impl( const state_type &x )
//	{
//	}
//
//	size_adjuster_type m_size_adjuster;
};




} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_EXPLICIT_EULER_HPP_INCLUDED
