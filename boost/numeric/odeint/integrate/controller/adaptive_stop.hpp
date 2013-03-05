/*
  [auto_generated]
  boost/numeric/odeint/integrate/controller/adaptive_stop.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_CONTROLLER_ADAPTIVE_STOP_HPP_DEFINED
#define BOOST_NUMERIC_ODEINT_INTEGRATE_CONTROLLER_ADAPTIVE_STOP_HPP_DEFINED

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>

#include <stdexcept>


namespace boost {
namespace numeric {
namespace odeint {


template< class Time = double >
class adaptive_stop
{
public:

    typedef Time time_type;

    adaptive_stop( time_type t_end ) : m_t_end( t_end ) { }

    template< class Stepper , class Sys , class State >
    void init( Stepper &stepper , Sys sys , const State &s , time_type t , time_type dt )
    {
    }

    template< class State >
    bool stop( const State &s , time_type &t )
    {
        return ( t > m_t_end ); 
    }

    template< class Stepper , class Sys , class State >
    void do_step( Stepper &stepper , Sys sys , State &x , time_type &t , time_type &dt )
    {
        typedef typename Stepper::stepper_category stepper_category;
        do_step_impl( stepper , sys , x , t , dt , stepper_category() );
    }

    template< class Stepper , class Sys , class State >
    void exit( Stepper &stepper , Sys sys , State &x , time_type t , time_type dt )
    {
    }

    



    template< class Stepper , class Sys , class State >
    void do_step_impl( Stepper &stepper , Sys sys , State &x , time_type &t , time_type &dt , stepper_tag )
    {
        stepper.do_step( sys , x , t , dt );
        t += dt;
    }

    template< class Stepper , class Sys , class State >
    void do_step_impl( Stepper &stepper , Sys sys , State &x , time_type &t , time_type &dt , controlled_stepper_tag )
    {
        size_t max_attempts = 1000;
        size_t trials = 0;
        controlled_step_result res = success;
        do
        {
            res = stepper.try_step( sys , x , t , dt );
            ++trials;
        }
        while( ( res == fail ) && ( trials < max_attempts ) );
        if( trials == max_attempts ) throw std::overflow_error( "adaptive_controller::do_step_impl : too much iterations" );

    }

    template< class Stepper , class Sys , class State >
    void do_step_impl( Stepper &stepper , Sys sys , State &x , time_type &t , time_type &dt , dense_output_stepper_tag )
    {
        stepper.do_step( sys );
        t = stepper.current_time();
    }



    bool m_stop;
    Time m_t_end;
};


} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_INTEGRATE_CONTROLLER_ADAPTIVE_STOP_HPP_DEFINED
