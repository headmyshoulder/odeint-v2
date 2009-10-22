/* Boost odeint/euler.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner
 
 This file includes the explicit euler solver for ordinary differential equations.

 It solves any ODE dx/dt = f(x,t) via
 x(t+dt) = x(t) + dt*f(x,t)

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_EULER_HPP
#define BOOST_NUMERIC_ODEINT_EULER_HPP

#include <boost/concept_check.hpp>

namespace boost {
namespace numeric {
namespace odeint {

    template< class ContainerType >
    class ode_step_euler
    {
	BOOST_CLASS_REQUIRE( ContainerType , boost , SequenceConcept );
	ContainerType dxdt;

    public:

	template< class DynamicalSystem , class TimeType>
	void next_step( DynamicalSystem system , ContainerType &x , TimeType t , TimeType dt )
	{
	    if( x.size() != dxdt.size() ) dxdt.resize( x.size() );
	    system( x , dxdt , t );
	    typename ContainerType::iterator state_begin = x.begin();
	    typename ContainerType::iterator state_end = x.end();
	    typename ContainerType::iterator derivative_begin = dxdt.begin();
	    while( state_begin != state_end ) 
		(*state_begin++) += dt * (*derivative_begin++);
	}
    };



/* ToDo:
   Write stepper for
   * fixed size systems
   * array<T>
   * system( T* , T* , T )
*/


} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_EULER_HPP
