/*
  [auto_generated]
  boost/numeric/odeint/iterator/detail/adaptive_time_iterator_dense_output_impl.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_ITERATOR_DETAIL_ADAPTIVE_TIME_ITERATOR_DENSE_OUTPUT_IMPL_HPP_DEFINED
#define BOOST_NUMERIC_ODEINT_ITERATOR_DETAIL_ADAPTIVE_TIME_ITERATOR_DENSE_OUTPUT_IMPL_HPP_DEFINED


namespace boost {
namespace numeric {
namespace odeint {

    /*
     * Specilization for steppers and error steppers
     */
    /**
     * \brief ODE Iterator with adaptive step size control. The value type of this iterator is a pair of state type and time type of the stepper.
     *
     * Implements an ODE with adaptive step size control. Uses dense-output steppers. adaptive_iterator is a model
     * of single-pass iterator.
     *
     * The value type of this iterator is a pair of state type and time type of the stepper.
     *
     * \tparam Stepper The stepper type which should be used during the iteration.
     * \tparam System The type of the system function (ODE) which should be solved.
     */
    template< class Stepper , class System , class State >
    class adaptive_time_iterator< Stepper , System , State , dense_output_stepper_tag > : public detail::ode_time_iterator_base
    <
        adaptive_time_iterator< Stepper , System , State , dense_output_stepper_tag > ,
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
            adaptive_time_iterator< Stepper , System , State , dense_output_stepper_tag > ,
            Stepper , System , State > base_type;


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
            : base_type( stepper , sys , s , t , dt ) , m_t_end( t_end )
        {
            if( detail::less_eq_with_sign( this->m_t , this->m_t_end , this->m_dt ) )
            {
                unwrapped_stepper_type &st = this->m_stepper;
                st.initialize( *( this->m_state ) , this->m_t , this->m_dt );
            } else {
                this->m_at_end = true;
            }
        }

        /**
         * \brief Constructs an adaptive_time_iterator. This constructor should be used to construct the end iterator.
         *
         * \param stepper The stepper to use during the iteration.
         * \param sys The system function (ODE) to solve.
         * \param s The initial state. adaptive_time_iterator stores a reference of s and changes its value during the iteration.
         */
        adaptive_time_iterator( stepper_type stepper , system_type sys , state_type &s )
            : base_type( stepper , sys , s )
        {
        }


    private:

        friend class boost::iterator_core_access;

        std::pair< const state_type& , const time_type& > dereference() const
        {
            const unwrapped_stepper_type &stepper = this->m_stepper;
            return std::pair< const state_type & , const time_type & >( stepper.current_state() , this->m_t );
        }


        void increment()
        {
            unwrapped_stepper_type &stepper = this->m_stepper;
            if( detail::less_with_sign( this->m_t ,
                                        this->m_t_end ,
                                        stepper.current_time_step() ) )
            {
                if( detail::less_with_sign( this->m_t_end ,
                                            static_cast<time_type>(this->m_t + stepper.current_time_step()) ,
                                            stepper.current_time_step() ) )
                {
                    // make stpper to end exactly at t_end
                    stepper.initialize( stepper.current_state() , stepper.current_time() ,
                                        static_cast<time_type>(this->m_t_end-this->m_t) );
                }
                stepper.do_step( this->m_system );
                this->m_t = stepper.current_time();
            } else { // we have reached t_end
                this->m_at_end = true;
            }
        }

    private:
        time_type m_t_end;
    };




} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_ITERATOR_DETAIL_ADAPTIVE_TIME_ITERATOR_DENSE_OUTPUT_IMPL_HPP_DEFINED
