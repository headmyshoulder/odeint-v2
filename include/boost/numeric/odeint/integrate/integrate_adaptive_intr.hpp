/*
 [auto_generated]
 boost/numeric/odeint/integrate/integrate_adaptive.hpp

 [begin_description]
 Adaptive integration of ODEs.
 [end_description]

 Copyright 2011-2013 Karsten Ahnert
 Copyright 2011-2012 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_ADAPTIVE_INTR_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_ADAPTIVE_INTR_HPP_INCLUDED

#include <boost/type_traits/is_same.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_adaptive_intr.hpp>

namespace boost {
namespace numeric {
namespace odeint {


/*
 * the two overloads are needed in order to solve the forwarding problem
 */
template< class Stepper , class System , class State , class Time , class Interrupt, class Observer >
size_t integrate_adaptive_interruptible(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Time end_time , Time dt , Interrupt interrupt,
        Observer observer )
{
    typedef typename odeint::unwrap_reference< Stepper >::type::stepper_category stepper_category;

    return detail::integrate_adaptive_interruptible(
            stepper , system , start_state ,
            start_time , end_time , dt , interrupt,
            observer , stepper_category() );
}

/************* DOXYGEN ************/

    /** 
     * \fn integrate_adaptive_interruptible( Stepper stepper , System system , State &start_state , Time start_time , Time end_time , Time dt , Interrupt interrupt , Observer observer )
     * \brief Integrates the ODE with adaptive step size.
     * 
     * This function integrates the ODE given by system with the given stepper.
     * This function is similar to integrate_adaptive except for allowing user to 
     * pass function interrupt which force some additional points into the solver.
     *
     * \param stepper The stepper to be used for numerical integration.
     * \param system Function/Functor defining the rhs of the ODE.
     * \param start_state The initial condition x0.
     * \param start_time The initial time t0.
     * \param end_time The final integration time tend.
     * \param dt The time step between observer calls, _not_ necessarily the 
     * time step of the integration.
     * \param interrupt Function/Functor called at forced point and its return value determines the next forced point.
     * \param observer Function/Functor called at equidistant time intervals.
     * \return The number of steps performed.
     */

} // namespace odeint
} // namespace numeric
} // namespace boost



#endif // BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_ADAPTIVE_HPP_INCLUDED
