/*
 boost header: NUMERIC_ODEINT_STEPPER/controlled_error_stepper.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_ERROR_STEPPER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_ERROR_STEPPER_HPP_INCLUDED

#include <boost/numeric/odeint/algebra/standard_algebra.hpp>
#include <boost/numeric/odeint/algebra/standard_operations.hpp>

#include <boost/numeric/odeint/stepper/explicit_stepper_base.hpp>
#include <boost/numeric/odeint/stepper/detail/macros.hpp>

namespace boost {
namespace numeric {
namespace odeint {


//template<
//    class State ,
//    class Time = double ,
//	class Algebra = standard_algebra< State > ,
//	class Operations = standard_operations< Time > ,
//	class AdjustSizePolicy = adjust_size_initially_tag
//	>
//class explicit_euler
//: public explicit_stepper_base<
//	  explicit_euler< State , Time , Algebra , Operations , AdjustSizePolicy > ,
//	  1 , State , Time , Algebra , Operations , AdjustSizePolicy >
//{
//public :
//
//	BOOST_ODEINT_EXPLICIT_STEPPERS_TYPEDEFS( explicit_euler , 1 );
//
//	template< class System >
//	void do_step_impl( System system , state_type &x , const state_type &dxdt , time_type t , time_type dt )
//	{
//		algebra_type::for_each2( x , dxdt , typename operations_type::increment1( dt ) );
//	}
//};




} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_ERROR_STEPPER_HPP_INCLUDED
