/*
  [auto_generated]
  boost/numeric/odeint/integrate/controller/approximate.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_CONTROLLER_APPROXIMATE_HPP_DEFINED
#define BOOST_NUMERIC_ODEINT_INTEGRATE_CONTROLLER_APPROXIMATE_HPP_DEFINED

#include <boost/numeric/odeint/util/copy.hpp>
#include <boost/numeric/odeint/util/resize.hpp>



namespace boost {
namespace numeric {
namespace odeint {



template< class Pred , class Time = double >
class approximate
{
public:


    approximate( Pred pred , Time min_step )
        : m_pred( pred ) , m_min_step( min_step ) , m_stop( false ) { }

    template< class Stepper , class Sys , class State >
    void init( Stepper &stepper , Sys sys , const State &s , Time t , Time dt )
    {
    }

    template< class State >
    bool stop( const State &s , Time &t )
    {
        return m_stop;
    }

    template< class Stepper , class Sys , class State >
    void do_step( Stepper &stepper , Sys sys , State &x , Time &t , Time &dt )
    {
        typename Stepper::state_type x_old( x );
        Time t_old = t;
        typedef typename Stepper::stepper_category stepper_category;
        do_step_impl( stepper , sys , x , t , dt , stepper_category() );
        if( m_pred( x , t ) )
        {
            t = t_old;
            boost::numeric::odeint::copy( x_old , x );
            dt *= 0.5;
        }
        if( dt < m_min_step ) m_stop = true;
    }

    template< class Stepper , class Sys , class State >
    void exit( Stepper &stepper , Sys sys , State &x , Time t , Time dt )
    {
    }

    



    template< class Stepper , class Sys , class State >
    void do_step_impl( Stepper &stepper , Sys sys , State &x , Time &t , Time &dt , stepper_tag )
    {
        stepper.do_step( sys , x , t , dt );
        t += dt;
    }

    template< class Stepper , class Sys , class State >
    void do_step_impl( Stepper &stepper , Sys sys , State &x , Time &t , Time &dt , controlled_stepper_tag )
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


    Pred m_pred ;
    Time m_min_step ;
    bool m_stop;
};


template< class Pred , class Time >
approximate< Pred , Time > make_approximate( Pred pred , Time min_step )
{
    return approximate< Pred , Time >( pred , min_step );
}



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_INTEGRATE_CONTROLLER_APPROXIMATE_HPP_DEFINED
