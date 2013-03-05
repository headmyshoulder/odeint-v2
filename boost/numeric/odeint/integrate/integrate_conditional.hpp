/*
  [auto_generated]
  boost/numeric/odeint/integrate/integrate_conditional.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_CONDITIONAL_HPP_DEFINED
#define BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_CONDITIONAL_HPP_DEFINED


namespace boost {
namespace numeric {
namespace odeint {


template< class Stepper , class Sys , class State , class Time , class Controller , class Observer >
void integrate_conditional( Stepper stepper , Sys sys , State &x , Time t , Time dt , Controller controller , Observer observer )
{
    controller.init( stepper , sys , x , t , dt );
    while( !controller.stop( x , t ) )
    {
        observer( x , t );
        controller.do_step( stepper , sys , x , t , dt );
    }
    controller.exit( stepper , sys , x , t , dt );
}





} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_CONDITIONAL_HPP_DEFINED
