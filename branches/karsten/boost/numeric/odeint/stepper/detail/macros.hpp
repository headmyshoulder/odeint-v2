/*
 boost header: BOOST_NUMERIC_ODEINT_DETAIL/macros.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_DETAIL_MACROS_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_DETAIL_MACROS_HPP_INCLUDED


#define BOOST_ODEINT_EXPLICIT_STEPPERS_TYPEDEFS( STEPPER , ORDER ) \
typedef explicit_stepper_base< \
STEPPER< State , Value , Deriv ,Time , Algebra , Operations , AdjustSizePolicy > , \
ORDER , State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > stepper_base_type; \
typedef typename stepper_base_type::state_type state_type; \
typedef typename stepper_base_type::value_type value_type; \
typedef typename stepper_base_type::deriv_type deriv_type; \
typedef typename stepper_base_type::time_type time_type; \
typedef typename stepper_base_type::algebra_type algebra_type; \
typedef typename stepper_base_type::operations_type operations_type; \
typedef typename stepper_base_type::adjust_size_policy adjust_size_policy; \
typedef typename stepper_base_type::stepper_type stepper_type


#define BOOST_ODEINT_EXPLICIT_ERROR_STEPPERS_TYPEDEFS( STEPPER , STEPPER_ORDER , ERROR_ORDER ) \
typedef explicit_error_stepper_base< \
STEPPER< State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > , \
STEPPER_ORDER , ERROR_ORDER , State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > stepper_base_type; \
typedef typename stepper_base_type::state_type state_type; \
typedef typename stepper_base_type::value_type value_type; \
typedef typename stepper_base_type::deriv_type deriv_type; \
typedef typename stepper_base_type::time_type time_type; \
typedef typename stepper_base_type::algebra_type algebra_type; \
typedef typename stepper_base_type::operations_type operations_type; \
typedef typename stepper_base_type::adjust_size_policy adjust_size_policy; \
typedef typename stepper_base_type::stepper_type stepper_type


#define BOOST_ODEINT_EXPLICIT_STEPPERS_AND_ERROR_STEPPERS_TYPEDEFS( STEPPER , ORDER , STEPPER_ORDER , ERROR_ORDER ) \
typedef explicit_stepper_and_error_stepper_base< \
STEPPER< State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > , \
ORDER , STEPPER_ORDER , ERROR_ORDER , State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > stepper_base_type; \
typedef typename stepper_base_type::state_type state_type; \
typedef typename stepper_base_type::value_type value_type; \
typedef typename stepper_base_type::deriv_type deriv_type; \
typedef typename stepper_base_type::time_type time_type; \
typedef typename stepper_base_type::algebra_type algebra_type; \
typedef typename stepper_base_type::operations_type operations_type; \
typedef typename stepper_base_type::adjust_size_policy adjust_size_policy; \
typedef typename stepper_base_type::stepper_type stepper_type

#define BOOST_ODEINT_EXPLICIT_STEPPERS_AND_ERROR_STEPPERS_FSAL_TYPEDEFS( STEPPER , ORDER , STEPPER_ORDER , ERROR_ORDER ) \
typedef explicit_stepper_and_error_stepper_fsal_base< \
STEPPER< State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > , \
ORDER , STEPPER_ORDER , ERROR_ORDER , State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > stepper_base_type; \
typedef typename stepper_base_type::state_type state_type; \
typedef typename stepper_base_type::value_type value_type; \
typedef typename stepper_base_type::deriv_type deriv_type; \
typedef typename stepper_base_type::time_type time_type; \
typedef typename stepper_base_type::algebra_type algebra_type; \
typedef typename stepper_base_type::operations_type operations_type; \
typedef typename stepper_base_type::adjust_size_policy adjust_size_policy; \
typedef typename stepper_base_type::stepper_type stepper_type


#endif //BOOST_BOOST_NUMERIC_ODEINT_DETAIL_MACROS_HPP_INCLUDED
