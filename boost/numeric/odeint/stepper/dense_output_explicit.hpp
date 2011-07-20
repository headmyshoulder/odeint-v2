/*
 [auto_generated]
 boost/numeric/odeint/stepper/dense_output_explicit.hpp

 [begin_description]
 Dense-output functionality for explicit steppers.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DENSE_OUTPUT_EXPLICIT_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_DENSE_OUTPUT_EXPLICIT_HPP_INCLUDED


#include <utility>

#include <boost/numeric/odeint/util/copy.hpp>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

namespace boost {
namespace numeric {
namespace odeint {


template
<
class Stepper
>
class dense_output_explicit
{

private:

    void copy_pointers( const dense_output_explicit &dense_output )
    {
        if( dense_output.m_current_state == (&dense_output.m_x1.m_v ) )
        {
            m_current_state = &m_x1.m_v;
            m_old_state = &m_x2.m_v;
        }
        else
        {
            m_current_state = &m_x2.m_v;
            m_old_state = &m_x1.m_v;
        }
    }

public:

    /*
     * We do not need all typedefs.
     */
    typedef Stepper stepper_type;
    typedef typename stepper_type::state_type state_type;
    typedef typename stepper_type::wrapped_state_type wrapped_state_type;
    typedef typename stepper_type::value_type value_type;
    typedef typename stepper_type::deriv_type deriv_type;
    typedef typename stepper_type::wrapped_deriv_type wrapped_deriv_type;
    typedef typename stepper_type::time_type time_type;
    typedef typename stepper_type::algebra_type algebra_type;
    typedef typename stepper_type::operations_type operations_type;
    typedef typename stepper_type::resizer_type resizer_type;
    typedef dense_output_stepper_tag stepper_category;
    typedef dense_output_explicit< Stepper > dense_output_stepper_type;



    dense_output_explicit( const stepper_type &stepper = stepper_type() )
    : m_stepper( stepper ) ,
      m_current_state( &m_x1.m_v ) , m_old_state( &m_x2.m_v )
    { }


    dense_output_explicit( const dense_output_explicit &dense_output )
    : m_stepper( dense_output.m_stepper ) , m_x1( dense_output.m_x1 ) , m_x2( dense_output.m_x2 ) ,
      m_current_state( &m_x1.m_v ) , m_old_state( &m_x2.m_v ) ,
      m_t( dense_output.m_t ) , m_t_old( dense_output.m_t_old ) , m_dt( dense_output.m_dt )
    {
        copy_pointers( dense_output );
    }

    dense_output_explicit& operator=( const dense_output_explicit &dense_output )
    {
        m_stepper = dense_output.m_stepper;
        m_x1 = dense_output.m_x1;
        m_x2 = dense_output.m_x2;
        m_t = dense_output.m_t;
        m_t_old = dense_output.m_t_old;
        m_dt = dense_output.m_dt;
        copy_pointers( dense_output );
        return *this;
    }

    template< class StateType >
    void initialize( const StateType &x0 , const time_type &t0 , const time_type &dt0 )
    {
        m_resizer.adjust_size( x0 , boost::bind( &dense_output_stepper_type::template resize< StateType > , boost::ref( *this ) , _1 ) );
        boost::numeric::odeint::copy( x0 , *m_current_state );
        m_t = t0;
        m_dt = dt0;
    }

    template< class System >
    std::pair< time_type , time_type > do_step( System system )
    {
        m_stepper.do_step( system , *m_current_state , m_t , *m_old_state , m_dt );
        m_t_old = m_t;
        m_t += m_dt;
        std::swap( m_current_state , m_old_state );
        return std::make_pair( m_t_old , m_dt );
    }

    /*
     * The next two overloads are needed to solve the forwarding problem
     */
    template< class StateOut >
    void calc_state( const time_type &t , StateOut &x )
    {
        m_stepper.calc_state( x , t , *m_old_state , m_t_old , *m_current_state , m_t );
    }

    template< class StateOut >
    void calc_state( const time_type &t , const StateOut &x )
    {
        m_stepper.calc_state( x , t , *m_old_state , m_t_old , *m_current_state , m_t );
    }

    template< class StateIn >
    bool resize( const StateIn &x )
    {
        bool resized = false;
        resized |= adjust_size_by_resizeability( m_x1 , x , typename wrapped_state_type::is_resizeable() );
        resized |= adjust_size_by_resizeability( m_x2 , x , typename wrapped_state_type::is_resizeable() );
        return resized;
    }

    template< class StateType >
    void adjust_size( const StateType &x )
    {
        resize( x );
        m_stepper.stepper().resize( x );
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


private:

    stepper_type m_stepper;
    resizer_type m_resizer;
    wrapped_state_type m_x1 , m_x2;
    state_type *m_current_state , *m_old_state;
    time_type m_t , m_t_old , m_dt;

};

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_ODEINT_STEPPER_DENSE_OUTPUT_EXPLICIT_HPP_INCLUDED
