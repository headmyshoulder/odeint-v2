/*
  [auto_generated]
  boost/numeric/odeint/integrate/controller/conditional_stop.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_CONTROLLER_CONDITIONAL_STOP_HPP_DEFINED
#define BOOST_NUMERIC_ODEINT_INTEGRATE_CONTROLLER_CONDITIONAL_STOP_HPP_DEFINED

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>

#include <stdexcept>


namespace boost {
namespace numeric {
namespace odeint {


template< class Pred >
class conditional_stop
{
public:

    conditional_stop( Pred pred ) : m_pred( pred ) { }

    template< class Stepper , class Sys , class State , class Time >
    void init( Stepper &stepper , Sys sys , const State &s , Time t , Time dt ) const
    {
    }

    template< class State , class Time >
    bool stop( const State &s , Time &t ) const
    {
        return m_pred( s , t );
    }

    template< class Stepper , class Sys , class State , class Time >
    void do_step( Stepper &stepper , Sys sys , State &x , Time &t , Time &dt ) const
    {
        typedef typename Stepper::stepper_category stepper_category;
        do_step_impl( stepper , sys , x , t , dt , stepper_category() );
    }

    template< class Stepper , class Sys , class State , class Time >
    void exit( Stepper &stepper , Sys sys , State &x , Time t , Time dt ) const
    {
    }

    



    template< class Stepper , class Sys , class State , class Time >
    void do_step_impl( Stepper &stepper , Sys sys , State &x , Time &t , Time &dt , stepper_tag ) const
    {
        stepper.do_step( sys , x , t , dt );
        t += dt;
    }

    template< class Stepper , class Sys , class State , class Time >
    void do_step_impl( Stepper &stepper , Sys sys , State &x , Time &t , Time &dt , controlled_stepper_tag ) const
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
    

    Pred m_pred;
};


template< class Pred >
conditional_stop< Pred > make_conditional_stop( Pred pred )
{
    return conditional_stop< Pred >( pred );
}


} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_INTEGRATE_CONTROLLER_CONDITIONAL_STOP_HPP_DEFINED
