
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

#include <boost/numeric/odeint/util/unit_helper.hpp>


namespace boost {
namespace numeric {
namespace odeint {
namespace detail {




    template< class Iterator , class Stepper , class System , class StepperTag >
    class ode_time_iterator_base : public boost::iterator_facade
    <
        Iterator ,
        std::pair< typename Stepper::state_type& , typename Stepper::time_type > const ,
        boost::single_pass_traversal_tag
    >
    {
    private:

        typedef Stepper stepper_type;
        typedef System system_type;
        typedef typename stepper_type::state_type state_type;
        typedef typename stepper_type::time_type time_type;
        typedef typename stepper_type::value_type ode_value_type;

    public:
   
        ode_time_iterator_base( stepper_type stepper , system_type sys , state_type &s , time_type t , time_type t_end , time_type dt )
            : m_stepper( stepper ) , m_system( sys ) , m_state( s , t ) , m_t_end( t_end ) , m_dt( dt ) , m_first( true )
        {
            check_end();
        }

        ode_time_iterator_base( stepper_type stepper , system_type sys , state_type &s )
            : m_stepper( stepper ) , m_system( sys ) , m_state( s , 0.0 ) , m_t_end() , m_dt() , m_first( false )
        {
        }


    protected:

        friend class boost::iterator_core_access;

        bool equal( ode_time_iterator_base const& other ) const
        {
            if( m_first == other.m_first )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        const std::pair< state_type& , time_type >& dereference() const
        {
            return m_state;
        }

        void check_end( void )
        {
            if( get_unit_value( m_dt ) > static_cast< ode_value_type >( 0.0 ) )
            {
                if( m_state.second > m_t_end )
                {
                    m_first = false;
                }
            }
            else
            {
                if( m_state.second < m_t_end )
                {
                    m_first = false;
                }
            }
        }

        stepper_type m_stepper;
        system_type m_system;
        std::pair< state_type& , time_type > m_state;
        time_type m_t_end;
        time_type m_dt;
        bool m_first;

    };







} // namespace detail
} // namespace odeint
} // namespace numeric
} // namespace boost



#endif // BOOST_NUMERIC_ODEINT_ITERATOR_DETAIL_ODE_TIME_ITERATOR_BASE_HPP_INCLUDED
