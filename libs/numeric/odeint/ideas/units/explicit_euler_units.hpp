/*
 boost header: BOOST_NUMERIC_ODEINT/explicit_euler.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_EXPLICIT_EULER_UNITS_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_EXPLICIT_EULER_UNITS_HPP_INCLUDED

#include <boost/numeric/odeint/algebra/standard_algebra.hpp>

#include "explicit_stepper_base_units.hpp"
#include "standard_operations_units.hpp"

namespace boost {
namespace numeric {
namespace odeint {

template<
    class Deriv ,
    class Value = double ,
    class Time = double ,
	class Algebra = standard_algebra ,
	class Operations = standard_operations_units ,
	class AdjustSizePolicy = adjust_size_initially_tag
	>
class explicit_euler_units
: public explicit_stepper_base_units<
	  explicit_euler_units< Deriv , Value , Time , Algebra , Operations , AdjustSizePolicy > ,
	  1 , Deriv , Value , Time , Algebra , Operations , AdjustSizePolicy >
{
public :

	typedef explicit_stepper_base_units<
		explicit_euler_units< Deriv , Value , Time , Algebra , Operations , AdjustSizePolicy > ,
		1 , Deriv , Value , Time , Algebra , Operations , AdjustSizePolicy > stepper_base_type;
	typedef typename stepper_base_type::deriv_type deriv_type;
	typedef typename stepper_base_type::value_type value_type;
	typedef typename stepper_base_type::time_type time_type;
	typedef typename stepper_base_type::algebra_type algebra_type;
	typedef typename stepper_base_type::operations_type operations_type;
	typedef typename stepper_base_type::adjust_size_policy adjust_size_policy;
	typedef typename stepper_base_type::stepper_type stepper_type;


	template< class System , class State >
	void do_step_impl( System system , const State &in , const deriv_type &dxdt , const time_type t , State & out , const time_type dt )
	{
		typename algebra_type::for_each3()( out , in , dxdt , operations_type::make_scale_sum2( static_cast< value_type >( 1.0 ) , dt ) );
	}
};




} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_EXPLICIT_EULER_UNITS_HPP_INCLUDED
