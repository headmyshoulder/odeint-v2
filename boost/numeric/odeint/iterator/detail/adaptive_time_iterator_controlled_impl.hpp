/*
  [auto_generated]
  boost/numeric/odeint/iterator/detail/adaptive_time_iterator_controlled_impl.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_ITERATOR_DETAIL_ADAPTIVE_TIME_ITERATOR_CONTROLLED_IMPL_HPP_DEFINED
#define BOOST_NUMERIC_ODEINT_ITERATOR_DETAIL_ADAPTIVE_TIME_ITERATOR_CONTROLLED_IMPL_HPP_DEFINED

#include <boost/numeric/odeint/util/unit_helper.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/iterator/detail/ode_time_iterator_base.hpp>


namespace boost {
namespace numeric {
namespace odeint {


    /*
     * Specilization for controlled steppers
     */
    /**
     * \brief ODE Iterator with adaptive step size control. The value type of this iterator is a pair of state type and time type of the stepper.
     *
     * Implements an ODE iterator with adaptive step size control. Uses controlled steppers. adaptive_iterator is a model
     * of single-pass iterator.
     *
     * The value type of this iterator is a pair of state type and time type of the stepper.
     *
     * \tparam Stepper The stepper type which should be used during the iteration.
     * \tparam System The type of the system function (ODE) which should be solved.
     */
    template< class Stepper , class System >
    class adaptive_time_iterator< Stepper , System , controlled_stepper_tag > : public detail::ode_time_iterator_base
    <
        adaptive_time_iterator< Stepper , System , controlled_stepper_tag > ,
        Stepper , System
    >
    {
    private:

        typedef Stepper stepper_type;
        typedef System system_type;
        typedef typename boost::numeric::odeint::unwrap_reference< stepper_type >::type unwrapped_stepper_type;
        typedef typename unwrapped_stepper_type::state_type state_type;
        typedef typename unwrapped_stepper_type::time_type time_type;
        typedef typename unwrapped_stepper_type::value_type ode_value_type;
        typedef detail::ode_time_iterator_base<
            adaptive_time_iterator< Stepper , System , controlled_stepper_tag > ,
            Stepper , System > base_type;


    public:

        /**
         * \brief Constructs an adaptive_time_iterator. This constructor should be used to construct the begin iterator.
         *
         * \param stepper The stepper to use during the iteration.
         * \param sys The system function (ODE) to solve.
         * \param s The initial state. adaptive_time_iterator stores a reference of s and changes its value during the iteration.
         * \param t The initial time.
         * \param t_end The end time, at which the iteration should stop.
         * \param dt The initial time step.
         */
        adaptive_time_iterator( stepper_type stepper , system_type sys , state_type &s , time_type t , time_type t_end , time_type dt )
            : base_type( stepper , sys , s , t , t_end , dt ) {}

        /**
         * \brief Constructs an adaptive_time_iterator. This constructor should be used to construct the end iterator.
         *
         * \param stepper The stepper to use during the iteration.
         * \param sys The system function (ODE) to solve.
         * \param s The initial state. adaptive_time_iterator stores a reference of s and changes its value during the iteration.
         */
        adaptive_time_iterator( stepper_type stepper , system_type sys , state_type &s )
            : base_type( stepper , sys , s ) {}

    private:

        friend class boost::iterator_core_access;

        void increment()
        {
            unwrapped_stepper_type &stepper = this->m_stepper;
            const size_t max_attempts = 1000;
            size_t trials = 0;
            controlled_step_result res = success;
            do
            {
                res = stepper.try_step( this->m_system , *this->m_state , this->m_t , this->m_dt );
                ++trials;
            }
            while( ( res == fail ) && ( trials < max_attempts ) );
            if( trials == max_attempts )
            {
                throw std::overflow_error( "Adaptive iterator : Maximal number of iterations reached. A step size could not be found." );
            }
            this->check_end();
        }
    };



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_ITERATOR_DETAIL_ADAPTIVE_TIME_ITERATOR_CONTROLLED_IMPL_HPP_DEFINED
