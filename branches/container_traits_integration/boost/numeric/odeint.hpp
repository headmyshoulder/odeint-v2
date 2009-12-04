/* Boost odeint.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner
 
 This file includes *all* headers needed for integration of ordinary differential equations.
 
 
 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_HPP
#define BOOST_NUMERIC_ODEINT_HPP

#include <boost/numeric/odeint/stepper_euler.hpp>
#include <boost/numeric/odeint/stepper_rk4.hpp>
#include <boost/numeric/odeint/stepper_rk4_classical.hpp>
#include <boost/numeric/odeint/stepper_rk5_ck.hpp>
#include <boost/numeric/odeint/stepper_rk_generic.hpp>
#include <boost/numeric/odeint/stepper_half_step.hpp>
#include <boost/numeric/odeint/stepper_midpoint.hpp>
#include <boost/numeric/odeint/stepper_rk78_fehlberg.hpp>

#include <boost/numeric/odeint/controlled_stepper_standard.hpp>
#include <boost/numeric/odeint/controlled_stepper_bs.hpp>

#include <boost/numeric/odeint/error_checker_standard.hpp>

#include <boost/numeric/odeint/integrator_adaptive_stepsize.hpp>
#include <boost/numeric/odeint/integrator_constant_stepsize.hpp>

#endif // BOOST_NUMERIC_ODEINT_HPP
