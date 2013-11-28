/*
  [auto_generated]
  boost/numeric/odeint/iterator/detail/const_step_time_iterator_impl.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_ITERATOR_DETAIL_CONST_STEP_TIME_ITERATOR_IMPL_HPP_DEFINED
#define BOOST_NUMERIC_ODEINT_ITERATOR_DETAIL_CONST_STEP_TIME_ITERATOR_IMPL_HPP_DEFINED

#include <boost/numeric/odeint/util/unit_helper.hpp>
#include <boost/numeric/odeint/iterator/detail/ode_time_iterator_base.hpp>


namespace boost {
namespace numeric {
namespace odeint {


    /*
     * Specilization for steppers and error steppers
     */
    /**
     * \brief ODE Iterator with constant step size. The value type of this iterator is a pair of state type and
     * time type of the stepper.
     *
     * Implements an ODE iterator which solves the ODE with constant step size. The stepper must fulfill the Stepper
     * concept. const_step_time_iterator is a model of single-pass iterator.
     *
     * The value type of this iterator is a pair of state type and time type of the stepper.
     *
     * \tparam Stepper The stepper type which should be used during the iteration.
     * \tparam System The type of the system function (ODE) which should be solved.
     */
    template< class Stepper , class System , class State >
    class const_step_time_iterator< Stepper , System , State , stepper_tag > : public detail::ode_time_iterator_base
    <
        const_step_time_iterator< Stepper , System , State , stepper_tag > ,
        Stepper , System , State
    >
    {
    private:

        typedef Stepper stepper_type;
        typedef System system_type;
        typedef typename boost::numeric::odeint::unwrap_reference< stepper_type >::type unwrapped_stepper_type;
        typedef State state_type;
        typedef typename unwrapped_stepper_type::time_type time_type;
        typedef typename unwrapped_stepper_type::value_type ode_value_type;
        typedef detail::ode_time_iterator_base<
            const_step_time_iterator< Stepper , System , State , stepper_tag > ,
            Stepper , System , State > base_type;

    public:
   
        /**
         * \brief Constructs a const_step_time_iterator. This constructor should be used to construct the begin iterator.
         *
         * \param stepper The stepper to use during the iteration.
         * \param sys The system function (ODE) to solve.
         * \param s The initial state. const_step_time_iterator stores a reference of s and changes its value during the iteration.
         * \param t The initial time.
         * \param t_end The end time, at which the iteration should stop.
         * \param dt The initial time step.
         */
        const_step_time_iterator( stepper_type stepper , system_type sys , state_type &s ,
                                  time_type t , time_type t_end , time_type dt )
            : base_type( stepper , sys , s , t , dt ) ,
              m_t_start( t ) , m_t_end( t_end ) , m_step(0)
        {
            if( detail::less_with_sign( this->m_t_end , this->m_t , this->m_dt ) )
                this->m_at_end = true;
        }

        /**
         * \brief Constructs a const_step_time_iterator. This constructor should be used to construct the end iterator.
         *
         * \param stepper The stepper to use during the iteration.
         * \param sys The system function (ODE) to solve.
         * \param s The initial state. const_step_time_iterator stores a reference of s and changes its value during the iteration.
         */
        const_step_time_iterator( stepper_type stepper , system_type sys , state_type &s )
            : base_type( stepper , sys , s ) , m_t_start() , m_t_end() , m_step()
        {}

    protected:

        friend class boost::iterator_core_access;

        void increment()
        {
            if( detail::less_eq_with_sign( static_cast<time_type>(this->m_t+this->m_dt) ,
                                           this->m_t_end , this->m_dt ) )
            {
                unwrapped_stepper_type &stepper = this->m_stepper;
                stepper.do_step( this->m_system , *this->m_state , this->m_t , this->m_dt );
                // use integer to compute current time to reduce roundoff errors
                this->m_step++;
                this->m_t = this->m_t_start + static_cast< typename unit_value_type<time_type>::type >(this->m_step)*this->m_dt;
            } else {
                this->m_at_end = true;
            }
        }

    private:
        time_type m_t_start;
        time_type m_t_end;
        size_t m_step;
    };

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_ITERATOR_DETAIL_CONST_STEP_TIME_ITERATOR_IMPL_HPP_DEFINED
