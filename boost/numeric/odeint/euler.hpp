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

#include <boost/numeric/odeint/detail/accumulators.hpp>
#include <boost/numeric/odeint/concepts/state_concept.hpp>
#include <boost/numeric/odeint/resizer.hpp>



namespace boost {
namespace numeric {
namespace odeint {

    template<
	class ContainerType ,
	class ResizeType = resizer< ContainerType >
	>
    class ode_step_euler
    {
	// check the concept of the ContainerType
        BOOST_CLASS_REQUIRE( ContainerType , boost::numeric::odeint, StateType );

        ContainerType dxdt;
	ContainerType xtemp;
        ResizeType resizer;

        typedef typename ContainerType::iterator iterator;
	typedef typename ContainerType::value_type value_type;

    public:

	template< class DynamicalSystem , class TimeType >
	void next_step( DynamicalSystem system ,
                        ContainerType &x ,
			ContainerType &dxdt ,
                        TimeType t ,
                        TimeType dt )
        {
	    detail::multiply_and_add( x.begin() , x.end() , dxdt.begin() , dt );
	}

        template< class DynamicalSystem , class TimeType >
        void next_step( DynamicalSystem system ,
                        ContainerType &x ,
                        TimeType t ,
                        TimeType dt )
        {
	    resizer.check_size_and_resize( x , dxdt );
            system( x , dxdt , t );
	    next_step( system , x , dxdt , t , dt );
        }


	template< class DynamicalSystem , class TimeType >
        void next_step( DynamicalSystem system ,
                        ContainerType &x ,
                        TimeType t ,
                        TimeType dt ,
			ContainerType &xerr )
        {
	    resizer.check_size_and_resize( x , dxdt );
	    resizer.check_size_and_resize( x , xerr );

	    xtemp = x;
	    TimeType dt2 = 0.5*dt;

	    system( x , dxdt , t );

	    next_step( system , x , dxdt , t , dt );
	    next_step( system , xtemp , dxdt , t , dt2 );
	    next_step( system , xtemp , t+dt2 , dt2 );

	    detail::substract_and_assign( x.begin() , x.end() , xtemp.begin() , xerr.begin() );
	}

    };

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_EULER_HPP
