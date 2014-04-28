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

#include <boost/iterator/iterator_facade.hpp>


namespace boost {
namespace numeric {
namespace odeint {
    
template< typename State , typename Stepper , typename Time >
class const_step_range
{
public:
    
    typedef State state_type;
    typedef Stepper stepper_type;
    typedef Time time_type;
    typedef const_step_range< state_type , stepper_type , time_type > range_type;
    
    struct const_step_iterator :
    public boost::iterator_facade
    <
        const_step_iterator ,
        state_type ,
        boost::single_pass_traversal_tag
    >
    {
        friend class boost::iterator_core_access;
        
        const_step_iterator( range_type &range , bool is_end )
        : m_range( range ) , m_is_end( is_end)
        { };
        
        
        void increment() { m_node = m_node->next(); }

        bool equal(node_iterator const& other) const
        {
            return this->m_node == other.m_node;
        }
        
        range_type &m_range;
        bool m_is_end;
    };
    
    typedef const_step_iterator iterator;
    
    const_step_range( state_type &x , stepper_type &stepper , time_type start_time , time_type end_time )
    : m_x( &x ) , m_stepper( &stepper ) , m_start_time( start_time ) , m_end_time( end_time )
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
    
    state_type *m_x;
    stepper_type *m_stepper;
    time_type m_start_time , m_end_time;
};

template< typename State , typename Stepper , typename Time >
const_step_range< State , Stepper , Time >
make_const_step_range( State &x , Stepper &stepper , Time start_time , Time end_time )
{
    return const_step_range< State , Stepper , Time >( x , stepper , start_time , end_time );
}




} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_RANGE_CONST_STEP_RANGE_HPP_INCLUDED
