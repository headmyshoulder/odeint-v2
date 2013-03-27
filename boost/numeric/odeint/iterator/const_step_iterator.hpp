
/*
 [auto_generated]
 boost/numeric/odeint/iterator/const_step_ode_iterator.hpp

 [begin_description]
 Iterator for iterating through the solution of an ODE with constant step size.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_ITERATOR_CONST_STEP_ODE_ITERATOR_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_ITERATOR_CONST_STEP_ODE_ITERATOR_HPP_INCLUDED


#include <boost/numeric/odeint/util/stepper_traits.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>



namespace boost {
namespace numeric {
namespace odeint {







    template< class Stepper , class System ,
              class StepperTag = typename base_tag< typename traits::stepper_category< Stepper >::type >::type >
    class const_step_iterator;






    template< class Stepper , class System >
    const_step_iterator< Stepper , System > make_const_step_iterator_begin(
        Stepper stepper ,
        System system , 
        typename traits::state_type< Stepper >::type &x ,
        typename traits::time_type< Stepper >::type t ,
        typename traits::time_type< Stepper >::type t_end ,
        typename traits::time_type< Stepper >::type dt )
    {
        return const_step_iterator< Stepper , System >( stepper , system , x , t , t_end , dt );
    }

    template< class Stepper , class System >
    const_step_iterator< Stepper , System > make_const_step_iterator_end(
        Stepper stepper ,
        System system , 
        typename traits::state_type< Stepper >::type &x )
    {
        return const_step_iterator< Stepper , System >( stepper , system , x );
    }

    template< class Stepper , class System >
    std::pair< const_step_iterator< Stepper , System > , const_step_iterator< Stepper , System > >
    make_const_step_range(
        Stepper stepper ,
        System system , 
        typename traits::state_type< Stepper >::type &x ,
        typename traits::time_type< Stepper >::type t_start ,
        typename traits::time_type< Stepper >::type t_end ,
        typename traits::time_type< Stepper >::type dt )
    {
        return std::make_pair(
            const_step_iterator< Stepper , System >( stepper , system , x , t_start , t_end , dt ) ,
            const_step_iterator< Stepper , System >( stepper , system , x )
            );
    }





    /**
     * \fn make_const_step_iterator_begin( Stepper stepper , System system , typename Stepper::state_type &x , typename Stepper::time_type t , typename Stepper::time_type t_end , typename Stepper::time_type dt )
     *
     * \brief Factory function for const_step_iterator. Constructs a begin iterator.
     *
     * \param stepper The stepper to use during the iteration.
     * \param system The system function (ODE) to solve.
     * \param x The initial state. const_step_iterator stores a reference of s and changes its value during the iteration.
     * \param t The initial time.
     * \param t_end The end time, at which the iteration should stop.
     * \param dt The initial time step.
     * \returns The const step iterator.
     */


    /**
     * \fn make_const_step_iterator_end( Stepper stepper , System system , typename Stepper::state_type &x )
     * \brief Factory function for const_step_iterator. Constructs a end iterator.
     *
     * \param stepper The stepper to use during the iteration.
     * \param system The system function (ODE) to solve.
     * \param x The initial state. const_step_iterator stores a reference of s and changes its value during the iteration.
     * \returns The const_step_iterator.
     */


    /**
     * \fn make_const_step_range( Stepper stepper , System system , typename Stepper::state_type &x , typename Stepper::time_type t_start , typename Stepper::time_type t_end , typename Stepper::time_type dt )
     *
     * \brief Factory function to construct a single pass range of const step iterators. A range is here a pair
     * of const_step_iterator.
     *
     * \param stepper The stepper to use during the iteration.
     * \param system The system function (ODE) to solve.
     * \param x The initial state. const_step_iterator store a reference of s and changes its value during the iteration.
     * \param t The initial time.
     * \param t_end The end time, at which the iteration should stop.
     * \param dt The initial time step.
     * \returns The const step range.
     */


} // namespace odeint
} // namespace numeric
} // namespace boost


#include <boost/numeric/odeint/iterator/impl/const_step_iterator_impl.hpp>
#include <boost/numeric/odeint/iterator/impl/const_step_iterator_dense_output_impl.hpp>



#endif // BOOST_NUMERIC_ODEINT_ITERATOR_CONST_STEP_ODE_ITERATOR_HPP_INCLUDED
