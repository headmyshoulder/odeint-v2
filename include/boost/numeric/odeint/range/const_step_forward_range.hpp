/*
  [auto_generated]
  boost/numeric/odeint/range/const_step_forward_range.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_RANGE_CONST_STEP_FORWARD_RANGE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_RANGE_CONST_STEP_FORWARD_RANGE_HPP_INCLUDED

#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <deque>

#define ODEINT_DEBUG_ITERATOR 1
#ifdef ODEINT_DEBUG_ITERATOR

#include <ostream>
#include <iostream>
using namespace std;

#endif // ODEINT_DEBUG_ITERATOR



namespace boost {
namespace numeric {
namespace odeint {
namespace range {


template< typename Stepper , typename State , typename System , typename Time >
class const_step_forward_range
{
public:

    typedef Stepper stepper_type;
    typedef typename boost::numeric::odeint::unwrap_reference< stepper_type >::type unwrapped_stepper_type;
    typedef State state_type;
    typedef System system_type;
    typedef typename boost::numeric::odeint::unwrap_reference< system_type >::type unwrapped_system_type;
    typedef Time time_type;
    typedef const_step_forward_range< stepper_type , state_type , system_type , time_type > range_type;
    

    struct const_step_forward_iterator : public boost::iterator_facade
    <
        const_step_forward_iterator ,
        state_type const ,
        boost::forward_traversal_tag
    >
    {
        friend class boost::iterator_core_access;

        const_step_forward_iterator( void )
        : m_range( nullptr ) {}
        
        const_step_forward_iterator( range_type const* range )
        : m_range( range )
        { };

        ~const_step_forward_iterator( void )
        {
            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "iter destruct: " << this << "\n";
            #endif
            if( m_range != nullptr )
                m_range->m_manager.destroy( this );
        }

        const_step_forward_iterator( const_step_forward_iterator const& iter )
        : m_range( iter.m_range )
        {
            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "iter copy: " << this << " from " << &iter << "\n";
            #endif
            if( m_range != nullptr )
                m_range->m_manager.create_copy( &iter , this );
        }

        const_step_forward_iterator( const_step_forward_iterator&& iter )
        : m_range( iter.m_range )
        {
            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "iter move: " << this << " from " << &iter << "\n";
            #endif

            if( m_range != nullptr )
                m_range->m_manager.move_iter( &iter , this );
        }

        const_step_forward_iterator& operator=( const_step_forward_iterator const &iter )
        {
            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "iter copy assignment: " << this << " from " << &iter << "\n";
            #endif

            m_range = iter.m_range;
            if( m_range != nullptr )
                m_range->m_manager.create_copy( &iter , this );
            return *this;
        }

        const_step_forward_iterator& operator=( const_step_forward_iterator &&iter )
        {
            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "iter move assignment: " << this << " from " << &iter << "\n";
            #endif

            m_range = iter.m_range;
            if( m_range != nullptr )
                m_range->m_manager.move_iter( &iter , this );
            return *this;
        }

        void increment( void )
        {
            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "iter increment: " << this << "\n";
            #endif
            if( m_range != nullptr )
            {
                m_range->m_manager.advance( this );
                if( !m_range->m_manager.valid_iterator( this ) )
                    m_range = nullptr;
            }
        }

        bool equal( const_step_forward_iterator const& other) const
        {
            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "iter equal: " << this << " == " << &other << "\n";
            cout << "  ranges: " << m_range << " == " << other.m_range << "\n";
            #endif
            return ( m_range == other.m_range );
        }
        
        const state_type& dereference( void ) const
        {
            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "iter derefence: " << this << "\n";
            #endif
            return m_range->m_manager.dereference( this );
        }
        
        range_type const* m_range;
    };
    
    typedef const_step_forward_iterator iterator;
    typedef const_step_forward_iterator const_iterator;
    
    const_step_forward_range( stepper_type stepper , state_type const& x , system_type system , time_type start_time , time_type end_time , time_type dt )
        : m_manager( this , x )
    , m_stepper( stepper )
    , m_system( system )
    , m_current_time( start_time ) , m_end_time( end_time ) , m_dt( dt )
    { }
    
    iterator begin( void )
    {
        return m_manager.new_iter();
    }

    iterator end( void )
    {
        return iterator( nullptr );
    }
    
    const_iterator begin( void ) const
    {
        return m_manager.new_iter();
    }
    
    const_iterator end( void ) const
    {
        return const_iterator( nullptr );
    }

    #ifdef ODEINT_DEBUG_ITERATOR
    void print_state( void ) const
    {
        m_manager.print_state();
    }
    #endif
    
    

    /**
     * HINTS: iterator needs move semantics - OK
     */
    struct iterator_manager
    {

        typedef std::vector< std::pair< iterator* , size_t > > positions_type;
        typedef std::deque< std::pair< size_t , state_type > > state_vector_type;
        typedef boost::shared_ptr< state_vector_type > state_vector_pointer;


        iterator_manager( range_type const* range , state_type const &x )
        : m_range( range )
        {
            m_states->push_back( std::make_pair( 0 , x ) );
        }

        #ifdef ODEINT_DEBUG_ITERATOR
        static std::string get_indent( size_t indent )
        {
            std::string str = "";
            for( size_t i=0 ; i<indent ; ++i )
                str += "  ";
            return str;
        }
        void print_state( size_t indent = 0 ) const
        {
            cout << get_indent( indent ) << "Manager state:" << "\n";
            cout << get_indent( indent + 1 ) << "Positions:" << "\n";
            for( auto const& p : m_positions )
                cout << get_indent( indent + 2 ) << p.first << " : " << p.second << "\n";
            cout << get_indent( indent + 1 ) << "States:" << "\n";
            for( auto const& s : *m_states )
                cout << get_indent( indent + 2 ) << s.first << " : " << ( & ( s.second ) ) << "\n";
            cout << get_indent( indent + 1 ) << "Min index: " << m_min_index << "\n";
            cout << get_indent( indent + 1 ) << "Max index: " << m_max_index << "\n";
        }
        #endif

        void destroy( iterator *i )
        {
            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "manager:destroy: " << i << "\n";
            #endif

            auto p = find_iterator_mutable( i );
            size_t pos = p->second;
            m_positions.erase( p );

            auto s = find_state_mutable( pos );
            assert( s != m_states->end() );
            if( s->first == m_min_index )
            {
                // TODO: remove states which are not referenced
                m_states->erase( s );
                ++m_min_index;
            }

            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "Destroying iterator with " << pos << "\n";
            print_state( 1 );
            #endif
        }

        void advance( iterator *i )
        {
            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "manager::advance: " << i << "\n";
            #endif

            auto p = find_iterator_mutable( i );
            size_t pos = p->second;
            if( pos == m_max_index )
            {
                auto s = find_state( pos );
                m_states->push_back( std::make_pair( pos + 1 , state_type {} ) );
                auto &ss = m_states->back();
                // TODO: resize newly created state
                m_range->m_stepper.do_step( m_range->m_system ,
                                            s->second ,
                                            m_range->m_current_time ,
                                            ss.second ,
                                            m_range->m_dt );
                m_range->m_current_time += m_range->m_dt;
                ++m_max_index;
                p->second = pos + 1;

                #ifdef ODEINT_DEBUG_ITERATOR
                cout << "  Creating new state" << "\n";
                print_state( 1 );
                #endif
            }
            else
            {
                #ifdef ODEINT_DEBUG_ITERATOR
                cout << "  Iterating existing state" << "\n";
                print_state();
                #endif

            }
        }

        typename positions_type::const_iterator find_iterator( iterator const* i ) const
        {
            return boost::find_if( m_positions , [&i]( std::pair< iterator* , size_t > const& p ) -> bool {
                    return p.first == i ; } );
        }

        typename positions_type::iterator find_iterator_mutable( iterator const* i )
        {
            return boost::find_if( m_positions , [&i]( std::pair< iterator* , size_t > const& p ) -> bool {
                    return p.first == i ; } );
        }


        typename state_vector_type::const_iterator find_state( size_t i ) const
        {
            return boost::find_if( *m_states , [&i]( std::pair< size_t , state_type > const &p ) -> bool {
                    return p.first == i; } );
        }

        typename state_vector_type::iterator find_state_mutable( size_t i ) const
        {
            return boost::find_if( *m_states , [&i]( std::pair< size_t , state_type > const &p ) -> bool {
                    return p.first == i; } );
        }


        // returns if iterator is still valid, or if it is already at the end of the range
        bool valid_iterator( iterator const* i ) const
        {
            return find_iterator( i ) != m_positions.end();
        }

        bool valid_state( size_t i ) const
        {
            return find_state( i ) != m_states->end();
        }

        state_type const& dereference( iterator const* i )
        {
            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "manager::dereference: " << i << "\n";
            #endif

            assert( valid_iterator( i ) );
            auto p = find_iterator( i );
            assert( valid_state( p->second ) );
            auto s = find_state( p->second );

            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "  position: " << p->second << ", " << s->first << ", " << ( & ( s->second ) ) << "\n";
            #endif

            return s->second;
        }

        iterator new_iter( void )
        {
            iterator iter { m_range };
            m_positions.push_back( std::make_pair( &iter , 0 ) );
            return iter;
        }

        void move_iter( iterator const* from , iterator *to )
        {
            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "manager::move_iter: " << to << " from " << from << "\n";
            #endif 

            auto p = find_iterator_mutable( from );
            size_t pos = p->second;
            m_positions.erase( p );
            m_positions.push_back( std::make_pair( to , pos ) );
        }

        void create_copy( iterator const* from , iterator *to )
        {
            #ifdef ODEINT_DEBUG_ITERATOR
            cout << "manager::copy_iter: " << to << " from " << from << "\n";
            #endif 

            auto p = find_iterator( from );
            m_positions.push_back( std::make_pair( to , p->second ) );
        }

        range_type const* m_range;
        positions_type m_positions;
        state_vector_pointer m_states = boost::make_shared< state_vector_type >();
        size_t m_min_index = 0;
        size_t m_max_index = 0;
    };
    
    
private:

    // copy construction ?
    mutable iterator_manager m_manager;


    mutable stepper_type m_stepper;
    mutable system_type m_system;
    mutable time_type m_current_time;
    time_type m_end_time;
    time_type m_dt;
};

template< typename Stepper , typename State , typename System , typename Time >
const_step_forward_range< Stepper , State , System , Time >
make_const_step_forward_range( Stepper stepper , State &x , System system , Time start_time , Time end_time , Time dt )
{
    return const_step_forward_range< Stepper , State , System , Time >( stepper , x , system , start_time , end_time , dt );
}




} // namespace range
} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_RANGE_CONST_STEP_FORWARD_RANGE_HPP_INCLUDED
