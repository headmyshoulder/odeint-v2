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

/*
class c_string_range
: public range_facade<c_string_range>
{
friend range_access;
char const * sz_;
char const & current() const { return *sz_; }
bool done() const { return *sz_ == '\0'; }
void next() { ++sz_; }
public:
c_string_range() = default;
explicit c_string_range(char const *sz) : sz_(sz)
{
assert(sz != nullptr);
}
};
*/

 // #include <boost/numeric/odeint/range/detail/infinite_range.hpp>

#include <range/v3/all.hpp>


namespace boost {
namespace numeric {
namespace odeint {

template< typename Stepper , typename System , typename State , typename Time >
class const_step_range : public ranges::range_facade< const_step_range< Stepper , System , State , Time > , true >
{
public:

    using stepper_type = Stepper;
    using system_type = System;
    using state_type = State;
    using time_type = Time;

    const_step_range( void ) = default;

    const_step_range( stepper_type stepper , system_type system , state_type state ,
                      time_type time , time_type dt )
        : m_stepper( std::move( stepper ) )
        , m_system( std::move( system ) )
        , m_state( std::move( state ) )
        , m_time( time )
        , m_dt( dt )
    {}

private:

    friend ranges::range_access;

    state_type const& current( void ) const
    {
        return m_state;
    }
    
    bool done( void ) const
    {
        return false;
    }
  
    void next( void )
    {
        m_stepper.do_step( m_system , m_state , m_time , m_dt );
        m_time += m_dt;
    }

    stepper_type m_stepper;
    system_type m_system;
    state_type m_state;
    time_type m_time;
    time_type m_dt;
};

template< typename Stepper , typename System , typename State , typename Time >
auto make_const_step_range( Stepper stepper , System system , State state , Time time , Time dt )
{
    return const_step_range< Stepper , System , State , Time >(
        std::move( stepper ) , std::move( system ) , std::move( state ) , time , dt );
}




} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_RANGE_CONST_STEP_RANGE_HPP_INCLUDED
