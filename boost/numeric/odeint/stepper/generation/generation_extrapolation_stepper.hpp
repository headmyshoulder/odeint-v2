/*
 [auto_generated]
 boost/numeric/odeint/stepper/generation/generation_extrapolation_stepper.hpp

 [begin_description]
 Enable the factory functions for the controller extrapolations stepper.
 [end_description]

 Copyright 2009-2013 Karsten Ahnert
 Copyright 2009-2013 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_GENERATION_GENERATION_EXTRAPOLATION_STEPPER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_GENERATION_GENERATION_EXTRAPOLATION_STEPPER_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/extrapolation_stepper.hpp>
#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
 
namespace boost {
namespace numeric {
namespace odeint {

template< unsigned short Order , class State , class Value , class Deriv , class Time , class Algebra , class Operations , class Resize >
struct get_controller< extrapolation_stepper< Order , State , Value , Deriv , Time , Algebra , Operations , Resize > >
{
    typedef extrapolation_stepper< Order , State , Value , Deriv , Time , Algebra , Operations , Resize > stepper_type;
    typedef controlled_runge_kutta< stepper_type > type;
};


} // odeint
} // numeric
} // boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_GENERATION_GENERATION_EXTRAPOLATION_STEPPER_HPP_INCLUDED
