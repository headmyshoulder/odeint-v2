/*
  [auto_generated]
  boost/numeric/odeint/range/const_step_range.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_RANGE_CONST_STEP_RANGE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_RANGE_CONST_STEP_RANGE_HPP_INCLUDED

#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>

#include <boost/iterator/iterator_facade.hpp>




namespace boost {
namespace numeric {
namespace odeint {
    
template< typename Stepper , typename State , typename System , typename Time >
class const_step_range
{
public:

    typedef Stepper stepper_type;    
    typedef State state_type;
    typedef System system_type;
    typedef Time time_type;
    typedef const_step_range< stepper_type , state_type , system_type , time_type > range_type;
    
    struct const_step_iterator :
    public boost::iterator_facade
    <
        const_step_iterator ,
        state_type ,
        boost::single_pass_traversal_tag
    >
    {
        friend class boost::iterator_core_access;
        
        const_step_iterator( range_type &range , bool at_end )
        : m_range( range ) , m_at_end( at_end)
        { };
        
        
        void increment( void )
        {
            if( detail::less_eq_with_sign( static_cast<time_type>(this->m_t+this->m_dt) ,
                                           this->m_t_end , this->m_dt ) )
            {
                m_range.m_stepper->do_step( *( m_range.m_system ) , *( m_range.m_x ) , m_range.m_current_time , m_range.m_dt );
                m_range.m_current_time += m_range.m_dt;
            } 
            else
            {
                m_at_end = true;
            }
        }

        bool equal( const_step_iterator const& other) const
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
        
        const state_type& dereference( void ) const
        {
            return *( m_range.m_x );
        }
        
        range_type &m_range;
        bool m_at_end;
    };
    
    typedef const_step_iterator iterator;
    
    const_step_range( stepper_type &stepper , state_type &x , system_type &system , time_type start_time , time_type end_time , time_type dt )
    : m_stepper( &stepper ) , m_x( &x ) , m_system( &system ) , m_current_time( start_time ) , m_end_time( end_time ) , m_dt( dt )
    { }
    
    const_step_iterator begin( void )
    {
        return const_step_iterator( *this , false );
    }

    const_step_iterator end( void )
    {
        return const_step_iterator( *this , true );
    }
    
    
    
private:

    stepper_type *m_stepper;    
    state_type *m_x;
    system_type *m_system;
    time_type m_current_time , m_end_time;
    time_type m_dt;
};

template< typename Stepper , typename State , typename System , typename Time >
const_step_range< Stepper , State , System , Time >
make_const_step_range( Stepper &stepper , State &x , System &system , Time start_time , Time end_time )
{
    return const_step_range< Stepper , State , System , Time >( stepper , x , system , start_time , end_time );
}




} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_RANGE_CONST_STEP_RANGE_HPP_INCLUDED
