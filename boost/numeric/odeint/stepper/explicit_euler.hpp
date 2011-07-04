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

#include <boost/numeric/odeint/stepper/base/explicit_stepper_base.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/stepper/detail/macros.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template< class State , class Value , class Deriv , class Time , class Algebra , class Operations , class AdjustSizePolicy >
class dense_output_explicit_euler;

template<
    class State ,
    class Value = double ,
    class Deriv = State ,
    class Time = Value ,
	class Algebra = range_algebra ,
	class Operations = default_operations ,
	class AdjustSizePolicy = adjust_size_initially_tag
	>
class explicit_euler
: public explicit_stepper_base<
	  explicit_euler< State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > ,
	  1 , State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy >
{
public :

	friend class dense_output_explicit_euler< State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy >;

	BOOST_ODEINT_EXPLICIT_STEPPERS_TYPEDEFS( explicit_euler , 1 );

	template< class System , class StateIn , class DerivIn , class StateOut >
	void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt , const time_type &t , StateOut &out , const time_type &dt )
	{
		algebra_type::for_each3( out , in , dxdt , typename operations_type::template scale_sum2< value_type , time_type >( 1.0 , dt ) );

	}

	template< class StateOut , class StateIn1 , class StateIn2 >
	void calc_state( StateOut &x , const time_type &t ,  const StateIn1 &old_state , const time_type &t_old , const StateIn2 &current_state , const time_type &t_new )
	{
		time_type delta = t - t_old;
		algebra_type::for_each3( x , old_state , stepper_base_type::m_dxdt , typename operations_type::template scale_sum2< value_type , time_type >( 1.0 , delta ) );
	}

};




} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_EXPLICIT_EULER_HPP_INCLUDED
