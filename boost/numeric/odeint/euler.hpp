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


namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

    template< class DynamicalSystem ,
	      class StateIterator ,
	      class DerivativeIterator ,
	      class TimeType >
    void euler_imp( DynamicalSystem system ,
		    StateIterator state_begin ,
		    StateIterator state_end ,
		    DerivativeIterator derivative_begin ,
		    TimeType time_step ,
		    TimeType time )
    {
	system( state_begin , derivative_begin , time );
	while( state_begin != state_end ) (*state_begin++) += time_step * (*derivative_begin++);
    }

} // namespace odeint
} // namespace numeric
} // namespace boost
} // namespace detail






namespace boost {
namespace numeric {
namespace odeint {

    template< class ContainerType >
    class ode_step_euler
    {

	ContainerType dxdt;

    public:

	template< class DynamicalSystem , class TimeType>
	void next_step( DynamicalSystem system , ContainerType &x , TimeType t , TimeType dt )
	{
	    if( x.size() != dxdt.size() ) dxdt.resize( x.size() );

	    typename ContainerType::iterator dxdt_iter = dxdt.begin();
	    typename ContainerType::iterator x_iter = x.begin();
	    typename ContainerType::iterator x_end = x.end();

	    detail::euler_imp( system , x_iter , x_end , dxdt_iter , dt , t );
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
