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

#include <boost/config.hpp>

#include <boost/numeric/odeint/stepper/explicit_euler.hpp>
#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>

#include <boost/numeric/odeint/stepper/explicit_error_rk54_ck.hpp>
#include <boost/numeric/odeint/stepper/explicit_error_dopri5.hpp>

#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>

#include <boost/numeric/odeint/stepper/dense_output_explicit.hpp>
#include <boost/numeric/odeint/stepper/dense_output_controlled_explicit_fsal.hpp>

/*
 * Including this algebra slows down the compilation time
 */
// #include <boost/numeric/odeint/algebra/fusion_algebra.hpp>

#include <boost/numeric/odeint/integrate/integrate.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>


#endif // BOOST_NUMERIC_ODEINT_HPP
