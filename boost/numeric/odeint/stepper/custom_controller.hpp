/*
  [auto_generated]
  boost/numeric/odeint/stepper/custom_controller.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_CUSTOM_CONTROLLER_HPP_DEFINED
#define BOOST_NUMERIC_ODEINT_STEPPER_CUSTOM_CONTROLLER_HPP_DEFINED

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/util/copy.hpp>

namespace boost {
namespace numeric {
namespace odeint {



    template<
        class Pred ,
        class StepPolicy ,
        class Stepper ,
        class Resizer = typename Stepper::resizer_type ,
        class StepperCategory = typename Stepper::stepper_category >
    class custom_controller;



    template< class Pred , class StepPolicy , class Stepper >
    custom_controller< Pred , StepPolicy , Stepper > make_custom_controller(
        Pred pred , StepPolicy step_policy , const Stepper &stepper )
    {
        return custom_controller< Pred , StepPolicy , Stepper >( pred , step_policy , stepper );
    }



    template< class Pred , class StepPolicy , class Stepper , class Resizer >
    class custom_controller< Pred , StepPolicy , Stepper , Resizer , stepper_tag >
    {
    public:

        typedef Pred predicate_type;
        typedef StepPolicy step_policy_type;
        typedef Stepper stepper_type;
        typedef typename stepper_type::time_type time_type;
        typedef typename stepper_type::state_type state_type;
        typedef typename stepper_type::deriv_type deriv_type;
        typedef typename stepper_type::value_type value_type;
        typedef controlled_stepper_tag stepper_category;
        typedef Resizer resizer_type;
        typedef typename stepper_type::wrapped_state_type wrapped_state_type;


        custom_controller(
            predicate_type pred = predicate_type() ,
            step_policy_type step_policy = step_policy_type() ,
            const stepper_type &stepper = stepper_type() )
            : m_pred( pred ) , m_step_policy( step_policy ) , m_stepper( stepper )
        {
        }


        /*
         * Version 1 : try_step( sys , x , t , dt )
         *
         * The overloads are needed to solve the forwarding problem
         */
        template< class System , class StateInOut >
        controlled_step_result try_step( System system , StateInOut &x , time_type &t , time_type &dt )
        {
            return try_step_v1( system , x , t, dt );
        }

        template< class System , class StateInOut >
        controlled_step_result try_step( System system , const StateInOut &x , time_type &t , time_type &dt )
        {
            return try_step_v1( system , x , t, dt );
        }




        /*
         * Version 2 : try_step( sys , in , t , out , dt )
         *
         * This version does not solve the forwarding problem, boost.range can not be used.
         */
        template< class System , class StateIn , class StateOut >
        controlled_step_result try_step( System system , const StateIn &in , time_type &t , StateOut &out , time_type &dt )
        {
            m_stepper.do_step( system , in , t , out , dt );
            bool res = m_pred( in , out , t , dt );
            if( res )
            {
                t += dt;
                return success;
            }
            else
            {
                dt = m_step_policy( in , out , t , dt );
                return fail;
            }
        }
       
    private:

        template< class System , class StateInOut >
        controlled_step_result try_step_v1( System system , StateInOut &x , time_type &t , time_type &dt )
        {
            m_xnew_resizer.adjust_size( x , detail::bind( &custom_controller::template resize_m_xnew_impl< StateInOut > , detail::ref( *this ) , detail::_1 ) );
            controlled_step_result res = try_step( system , x , t , m_xnew.m_v , dt );
            if( res == success )
            {
                boost::numeric::odeint::copy( m_xnew.m_v , x );
            }
            return res;
        }

        template< class StateIn >
        bool resize_m_xnew_impl( const StateIn &x )
        {
            return adjust_size_by_resizeability( m_xnew , x , typename is_resizeable<state_type>::type() );
        }

        predicate_type m_pred;
        step_policy_type m_step_policy;
        stepper_type m_stepper;

        resizer_type m_xnew_resizer;
        wrapped_state_type m_xnew;
    };



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_CUSTOM_CONTROLLER_HPP_DEFINED
