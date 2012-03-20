/*
 [auto_generated]
 boost/numeric/odeint/stepper/generic_dense_output_stepper.hpp

 [begin_description]
 Implementation of a generic dense output stepper
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_GENERIC_DENSE_OUTPUT_STEPPER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_GENERIC_DENSE_OUTPUT_STEPPER_HPP_INCLUDED


#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/generic_dense_output_stepper_definition.hpp>

#include <boost/bind.hpp>



namespace boost {
namespace numeric {
namespace odeint {


template< class ControlledStepper >
class generic_dense_output_stepper< ControlledStepper , controlled_stepper_tag >
{
public:

    void copy_variables( const generic_dense_output_stepper &stepper )
    {
        m_stepper = stepper.m_stepper;
        m_x1 = stepper.m_x1;
        m_x2 = stepper.m_x2;
        if( stepper.m_current_state == ( & ( stepper.m_x1 ) ) )
        {
            m_current_state = &m_x1;
            m_old_state = &m_x2;
        }
        else
        {
            m_current_state = &m_x2;
            m_old_state = &m_x1;
        }
        m_t = stepper.m_t;
        m_t_old = stepper.m_t_old;
        m_dt = stepper.m_dt;
    }

public:

    typedef ControlledStepper controlled_stepper_type;
    typedef typename controlled_stepper_type::stepper_type stepper_type;
    typedef typename stepper_type::value_type value_type;
    typedef typename stepper_type::state_type state_type;
    typedef typename stepper_type::time_type time_type;
    typedef typename stepper_type::deriv_type deriv_type;
    typedef typename stepper_type::resizer_type resizer_type;
    typedef state_wrapper< state_type > wrapped_state_type;
    typedef state_wrapper< deriv_type > wrapped_deriv_type;

    typedef dense_output_stepper_tag stepper_category;

    typedef generic_dense_output_stepper< ControlledStepper > dense_output_stepper_type;

    generic_dense_output_stepper( const controlled_stepper_type &stepper = controlled_stepper_type() )
    : m_stepper( stepper ) ,
      m_x1() , m_x2() , m_current_state( &m_x1.m_v ) , m_old_state( &m_x2.m_v ) ,
      m_t() , m_t_old() , m_dt()
    { }

    generic_dense_output_stepper( const generic_dense_output_stepper &rb )
    : m_current_state( &m_x1.m_v ) , m_old_state( &m_x2.m_v )
    { }

    generic_dense_output_stepper& operator=( const generic_dense_output_stepper &rb )
    {
        copy_variables( rb );
        return *this;
    }



    template< class StateType >
    void initialize( const StateType &x0 , const time_type &t0 , const time_type &dt0 )
    {
        m_resizer.adjust_size( x0 , boost::bind( &dense_output_stepper_type::template resize_impl< StateType > , boost::ref( *this ) , _1 ) );
        *m_current_state = x0;
        m_t = t0;
        m_dt = dt0;
    }

    template< class System >
    std::pair< time_type , time_type > do_step( System system )
    {
        const size_t max_count = 1000;

        controlled_step_result res = fail;
        m_t_old = m_t;
        size_t count = 0;
        do
        {
            res = m_stepper.try_step( system , *m_current_state , m_t , *m_old_state , m_dt );
            if( count++ == max_count )
                throw std::overflow_error( "generic_dense_output_stepper : too much iterations!");
        }
        while( res == fail );
        m_stepper.stepper().prepare_dense_output();
        std::swap( m_current_state , m_old_state );
        return std::make_pair( m_t_old , m_t );
    }


    /*
     * The two overloads are needed in order to solve the forwarding problem.
     */
    template< class StateOut >
    void calc_state( const time_type &t , StateOut &x )
    {
        m_stepper.stepper().calc_state( t , x , *m_old_state , m_t_old , *m_current_state , m_t );
    }

    template< class StateOut >
    void calc_state( const time_type &t , const StateOut &x )
    {
        m_stepper.stepper().calc_state( t , x , *m_old_state , m_t_old , *m_current_state , m_t );
    }


    template< class StateType >
    void adjust_size( const StateType &x )
    {
        m_stepper.adjust_size( x );
        resize_impl( x );
    }




    const state_type& current_state( void ) const
    {
        return *m_current_state;
    }

    const time_type& current_time( void ) const
    {
        return m_t;
    }

    const time_type& previous_state( void ) const
    {
        return *m_old_state;
    }

    const time_type& previous_time( void ) const
    {
        return m_t_old;
    }

    const time_type& current_time_step( void ) const
    {
        return m_dt;
    }




private:

    template< class StateIn >
    bool resize_impl( const StateIn &x )
    {
        bool resized = false;
        resized |= adjust_size_by_resizeability( m_x1 , x , typename is_resizeable<state_type>::type() );
        resized |= adjust_size_by_resizeability( m_x2 , x , typename is_resizeable<state_type>::type() );
        return resized;
    }


    controlled_stepper_type m_stepper;
    resizer_type m_resizer;
    wrapped_state_type m_x1 , m_x2;
    state_type *m_current_state , *m_old_state;
    time_type m_t , m_t_old , m_dt;


};




}
}
}



#endif // BOOST_NUMERIC_ODEINT_STEPPER_GENERIC_DENSE_OUTPUT_STEPPER_HPP_INCLUDED
