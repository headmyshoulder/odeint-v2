
/*
 [auto_generated]
 boost/numeric/odeint/iterator/adaptive_time_iterator.hpp

 [begin_description]
 Iterator for iterating throught the solution of an ODE with adaptive step size. The dereferenced types containes also the time.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_ITERATOR_ADAPTIVE_TIME_ITERATOR_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_ITERATOR_ADAPTIVE_TIME_ITERATOR_HPP_INCLUDED

# include <boost/iterator/iterator_facade.hpp>

#include <boost/numeric/odeint/util/unit_helper.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>


namespace boost {
namespace numeric {
namespace odeint {


    template< class Stepper , class System , class StepperTag = typename base_tag< typename Stepper::stepper_category >::type > 
    class adaptive_time_iterator;


    /*
     * Specilization for controlled steppers
     */
    template< class Stepper , class System >
    class adaptive_time_iterator< Stepper , System , controlled_stepper_tag > : public boost::iterator_facade
    <
        adaptive_time_iterator< Stepper , System , controlled_stepper_tag > ,
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
   
        adaptive_time_iterator( stepper_type stepper , system_type sys , state_type &s , time_type t , time_type dt , bool first  )
            : m_stepper( stepper ) , m_system( sys ) , m_state( s , t ) , m_dt( dt ) , m_first( first ) {}

    private:

        friend class boost::iterator_core_access;

        void increment()
        {
            const size_t max_attempts = 1000;
            size_t trials = 0;
            controlled_step_result res = success;
            do
            {
                res = m_stepper.try_step( m_system , m_state.first , m_state.second , m_dt );
                ++trials;
            }
            while( ( res == fail ) && ( trials < max_attempts ) );
            if( trials == max_attempts )
            {
                throw std::overflow_error( "Adaptive iterator : Maximal number of iterations reached. A step size could not be found." );
            }
        }

        bool equal( adaptive_time_iterator const& other ) const
        {
            if( m_first == other.m_first )
            {
                return true;
            }
            else
            {
                if( m_first )
                {
                    return ( get_unit_value( m_dt ) > static_cast< ode_value_type >( 0.0 ) ) ?
                        ( m_state.second > other.m_state.second ) :
                        ( m_state.second < other.m_state.second ) ;
                }
                else
                {
                    return ( get_unit_value( m_dt ) > static_cast< ode_value_type >( 0.0 ) ) ?
                        ( m_state.second < other.m_state.second ) :
                        ( m_state.second > other.m_state.second ) ;
                }
            }
        }

        const std::pair< state_type& , time_type >& dereference() const
        {
            return m_state;
        }

        stepper_type m_stepper;
        system_type m_system;
        std::pair< state_type& , time_type > m_state;
        time_type m_dt;
        bool m_first;
    };






    /*
     * Specilization for steppers and error steppers
     */
    template< class Stepper , class System >
    class adaptive_time_iterator< Stepper , System , dense_output_stepper_tag > : public boost::iterator_facade
    <
        adaptive_time_iterator< Stepper , System , dense_output_stepper_tag > ,
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
   
        adaptive_time_iterator( stepper_type stepper , system_type sys , state_type &s , time_type t , time_type dt , bool first )
            : m_stepper( stepper ) , m_system( sys ) , m_state( s , t ) , m_dt( dt ) , m_first( first )
        {
            m_stepper.initialize( m_state.first , m_state.second , m_dt );
        }

    private:

        friend class boost::iterator_core_access;

        void increment()
        {
            m_state.second += m_dt;
            while(  m_stepper.current_time() < m_state.second )
                m_stepper.do_step( m_system );
            m_stepper.calc_state( m_state.second , m_state.first );
        }

        bool equal( adaptive_time_iterator const& other ) const
        {
            if( m_first == other.m_first )
            {
                return true;
            }
            else
            {
                if( m_first )
                {
                    return ( get_unit_value( m_dt ) > static_cast< ode_value_type >( 0.0 ) ) ?
                        ( m_state.second > other.m_state.second ) :
                        ( m_state.second < other.m_state.second ) ;
                }
                else
                {
                    return ( get_unit_value( m_dt ) > static_cast< ode_value_type >( 0.0 ) ) ?
                        ( m_state.second < other.m_state.second ) :
                        ( m_state.second > other.m_state.second ) ;
                }
            }
        }

        const std::pair< state_type& , time_type >& dereference() const
        {
            return m_state;
        }

        stepper_type m_stepper;
        system_type m_system;
        std::pair< state_type& , time_type > m_state;
        time_type m_dt;
        bool m_first;
    };








    template< class Stepper , class System >
    adaptive_time_iterator< Stepper , System > make_adaptive_time_iterator_begin(
        Stepper stepper ,
        System system , 
        typename Stepper::state_type &x ,
        typename Stepper::time_type t ,
        typename Stepper::time_type dt )
    {
        return adaptive_time_iterator< Stepper , System >( stepper , system , x , t , dt , true );
    }

    template< class Stepper , class System >
    adaptive_time_iterator< Stepper , System > make_adaptive_time_iterator_end(
        Stepper stepper ,
        System system , 
        typename Stepper::state_type &x ,
        typename Stepper::time_type t ,
        typename Stepper::time_type dt )
    {
        return adaptive_time_iterator< Stepper , System >( stepper , system , x , t , dt , false );
    }


    template< class Stepper , class System >
    std::pair< adaptive_time_iterator< Stepper , System > , adaptive_time_iterator< Stepper , System > >
    make_adaptive_time_range(
        Stepper stepper ,
        System system , 
        typename Stepper::state_type &x ,
        typename Stepper::time_type t_start ,
        typename Stepper::time_type t_end ,
        typename Stepper::time_type dt )
    {
        return std::make_pair(
            adaptive_time_iterator< Stepper , System >( stepper , system , x , t_start , dt , true ) ,
            adaptive_time_iterator< Stepper , System >( stepper , system , x , t_end , dt , false ) );
    }












} // namespace odeint
} // namespace numeric
} // namespace boost



#endif // BOOST_NUMERIC_ODEINT_ITERATOR_ADAPTIVE_TIME_ITERATOR_HPP_INCLUDED
