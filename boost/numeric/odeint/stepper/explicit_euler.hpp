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

#include <boost/numeric/odeint/stepper/base/explicit_stepper_base.hpp>
#include <boost/numeric/odeint/stepper/detail/macros.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template< class State , class Time , class Algebra , class Operations , class AdjustSizePolicy >
class dense_output_explicit_euler;

template<
    class State ,
    class Value = double ,
    class Deriv = State ,
    class Time = Value ,
	class Algebra = standard_algebra ,
	class Operations = standard_operations ,
	class AdjustSizePolicy = adjust_size_initially_tag
	>
class explicit_euler
: public explicit_stepper_base<
	  explicit_euler< State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > ,
	  1 , State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy >
{
public :

	friend class dense_output_explicit_euler< Deriv , Time , Algebra , Operations , AdjustSizePolicy >;

	BOOST_ODEINT_EXPLICIT_STEPPERS_TYPEDEFS( explicit_euler , 1 );

	template< class System , class StateIn , class DerivIn , class StateOut >
	void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt , const time_type &t , StateOut &out , const time_type &dt )
	{
		algebra_type::for_each3( out , in , dxdt , typename operations_type::template scale_sum2< value_type , time_type >( 1.0 , dt ) );
	}
};




} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_EXPLICIT_EULER_HPP_INCLUDED
