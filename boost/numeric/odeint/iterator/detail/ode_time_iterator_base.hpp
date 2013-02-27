
/*
 [auto_generated]
 boost/numeric/odeint/iterator/detail/ode_time_iterator_base.hpp

 [begin_description]
 Base class for const_step_time_iterator and adaptive_time_iterator.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_ITERATOR_DETAIL_ODE_TIME_ITERATOR_BASE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_ITERATOR_DETAIL_ODE_TIME_ITERATOR_BASE_HPP_INCLUDED

#include <boost/iterator/iterator_facade.hpp>

#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {




    template< class Iterator , class Stepper , class System >
    class ode_time_iterator_base : public boost::iterator_facade
    <
        Iterator ,
        std::pair< const typename traits::state_type< Stepper >::type , const typename traits::time_type< Stepper >::type > ,
        boost::single_pass_traversal_tag ,
        std::pair< const typename traits::state_type< Stepper >::type& , const typename traits::time_type< Stepper >::type& >
    >
    {
    private:

        typedef Stepper stepper_type;
        typedef System system_type;
        typedef typename boost::numeric::odeint::unwrap_reference< stepper_type >::type unwrapped_stepper_type;
        typedef typename unwrapped_stepper_type::state_type state_type;
        typedef typename unwrapped_stepper_type::time_type time_type;
        typedef typename unwrapped_stepper_type::value_type ode_value_type;

    public:
   
        ode_time_iterator_base( stepper_type stepper , system_type sys , state_type &s ,
                                time_type t , time_type t_end , time_type dt )
            : m_stepper( stepper ) , m_system( sys ) , m_state( &s ) , m_t( t ) , m_t_end( t_end ) , m_dt( dt ) , m_at_end( false )
        {
            check_end();
        }

        ode_time_iterator_base( stepper_type stepper , system_type sys , state_type &s )
            : m_stepper( stepper ) , m_system( sys ) , m_state( &s ) , m_t() , m_t_end() , m_dt() , m_at_end( true )
        {
        }

        bool same( ode_time_iterator_base const& iter )
        {
            return (
                ( m_state == iter.m_state ) &&
                ( m_t == iter.m_t ) &&
                ( m_t_end == iter.m_t_end ) &&
                ( m_dt == iter.m_dt ) &&
                ( m_at_end == iter.m_at_end )
                );
        }


    protected:

        friend class boost::iterator_core_access;

        bool equal( ode_time_iterator_base const& other ) const
        {
            if( m_at_end == other.m_at_end )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        std::pair< const state_type& , const time_type& > dereference() const
        {
            return std::pair< const state_type & , const time_type & >( *m_state , m_t );
        }

        void check_end( void )
        {
            if( detail::less_with_sign(  m_t_end , m_t , m_dt ) )
                m_at_end = true;
        }

        stepper_type m_stepper;
        system_type m_system;
        state_type *m_state;
        time_type m_t;
        time_type m_t_end;
        time_type m_dt;
        bool m_at_end;

    };







} // namespace detail
} // namespace odeint
} // namespace numeric
} // namespace boost



#endif // BOOST_NUMERIC_ODEINT_ITERATOR_DETAIL_ODE_TIME_ITERATOR_BASE_HPP_INCLUDED
