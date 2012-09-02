
/*
 [auto_generated]
 boost/numeric/odeint/iterator/const_step_ode_iterator.hpp

 [begin_description]
 Iterator for iterating throught the solution of an ODE with constant step size.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_ITERATOR_CONST_STEP_ODE_ITERATOR_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_ITERATOR_CONST_STEP_ODE_ITERATOR_HPP_INCLUDED

#include <boost/numeric/odeint/iterator/detail/ode_iterator_base.hpp>

#include <boost/numeric/odeint/util/unit_helper.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

namespace boost {
namespace numeric {
namespace odeint {


    template< class Stepper , class System , class StepperTag = typename base_tag< typename Stepper::stepper_category >::type > 
    class const_step_iterator;





    /*
     * Specilization for steppers and error steppers
     */
    template< class Stepper , class System >
    class const_step_iterator< Stepper , System , stepper_tag >
        : public detail::ode_iterator_base<
        const_step_iterator< Stepper , System , stepper_tag > ,
        Stepper , System , stepper_tag >
    {
    private:

        typedef Stepper stepper_type;
        typedef System system_type;
        typedef typename stepper_type::state_type state_type;
        typedef typename stepper_type::time_type time_type;
        typedef typename stepper_type::value_type ode_value_type;
        typedef detail::ode_iterator_base<
            const_step_iterator< Stepper , System , stepper_tag > ,
            Stepper , System , stepper_tag > base_type;

    public:
   
        const_step_iterator( stepper_type stepper , system_type sys , state_type &s , time_type t , time_type t_end , time_type dt )
            : base_type( stepper , sys , s , t , t_end , dt ) { }

        const_step_iterator( stepper_type stepper , system_type sys , state_type &s )
            : base_type( stepper , sys , s ) { }

    protected:

        friend class boost::iterator_core_access;

        void increment()
        {
            this->m_stepper.do_step( this->m_system , this->m_state , this->m_t , this->m_dt );
            this->m_t += this->m_dt;
            this->check_end();
        }
    };



    /*
     * Specilization for dense output stepper
     */
    template< class Stepper , class System >
    class const_step_iterator< Stepper , System , dense_output_stepper_tag >
        : public detail::ode_iterator_base<
        const_step_iterator< Stepper , System , dense_output_stepper_tag > ,
        Stepper , System , dense_output_stepper_tag >
    {
    private:

        typedef Stepper stepper_type;
        typedef System system_type;
        typedef typename stepper_type::state_type state_type;
        typedef typename stepper_type::time_type time_type;
        typedef typename stepper_type::value_type ode_value_type;
        typedef detail::ode_iterator_base<
            const_step_iterator< Stepper , System , dense_output_stepper_tag > ,
            Stepper , System , dense_output_stepper_tag > base_type;

    public:
   
        const_step_iterator( stepper_type stepper , system_type sys , state_type &s , time_type t , time_type t_end , time_type dt )
            : base_type( stepper , sys , s , t , t_end , dt )
        {
            this->m_stepper.initialize( this->m_state , this->m_t , this->m_dt );
        }

        const_step_iterator( stepper_type stepper , system_type sys , state_type &s )
            : base_type( stepper , sys , s ) { }


    protected:

        friend class boost::iterator_core_access;

        void increment( void )
        {
            this->m_t += this->m_dt;
            if( get_unit_value( this->m_dt ) > static_cast< ode_value_type >( 0.0 ) )
            {
                while( this->m_stepper.current_time() < this->m_t )
                    this->m_stepper.do_step( this->m_system );
            }
            else
            {
                while( this->m_stepper.current_time() > this->m_t )
                    this->m_stepper.do_step( this->m_system );
            }
            this->m_stepper.calc_state( this->m_t , this->m_state );
            this->check_end();
        }
    };





    template< class Stepper , class System >
    const_step_iterator< Stepper , System > make_const_step_iterator_begin(
        Stepper stepper ,
        System system , 
        typename Stepper::state_type &x ,
        typename Stepper::time_type t ,
        typename Stepper::time_type t_end ,
        typename Stepper::time_type dt )
    {
        return const_step_iterator< Stepper , System >( stepper , system , x , t , t_end , dt );
    }

    template< class Stepper , class System >
    const_step_iterator< Stepper , System > make_const_step_iterator_end(
        Stepper stepper ,
        System system , 
        typename Stepper::state_type &x )
    {
        return const_step_iterator< Stepper , System >( stepper , system , x );
    }


    template< class Stepper , class System >
    std::pair< const_step_iterator< Stepper , System > , const_step_iterator< Stepper , System > >
    make_const_step_range(
        Stepper stepper ,
        System system , 
        typename Stepper::state_type &x ,
        typename Stepper::time_type t_start ,
        typename Stepper::time_type t_end ,
        typename Stepper::time_type dt )
    {
        return std::make_pair(
            const_step_iterator< Stepper , System >( stepper , system , x , t_start , t_end , dt ) ,
            const_step_iterator< Stepper , System >( stepper , system , x )
            );
    }












} // namespace odeint
} // namespace numeric
} // namespace boost



#endif // BOOST_NUMERIC_ODEINT_ITERATOR_CONST_STEP_ODE_ITERATOR_HPP_INCLUDED
