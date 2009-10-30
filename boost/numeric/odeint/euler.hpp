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

#include <boost/numeric/odeint/detail/iterator_algebra.hpp>
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

	const unsigned int order() { return 1; }

	template< class DynamicalSystem , class TimeType >
	void next_step( DynamicalSystem system ,
                        ContainerType &x ,
			ContainerType &dxdt ,
                        TimeType t ,
                        TimeType dt )
        {
	    detail::it_algebra::scale_and_add( x.begin() , x.end() , dxdt.begin() , dt );
	}

        template< class DynamicalSystem , class TimeType >
        void next_step( DynamicalSystem system ,
                        ContainerType &x ,
                        TimeType t ,
                        TimeType dt )
        {
	    resizer.adjust_size( x , dxdt );
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
	    resizer.adjust_size( x , dxdt );
	    resizer.adjust_size( x , xerr );

	    xtemp = x;
	    TimeType dt2 = 0.5*dt;

	    system( x , dxdt , t );

	    next_step( system , x , dxdt , t , dt );
	    next_step( system , xtemp , dxdt , t , dt2 );
	    next_step( system , xtemp , t+dt2 , dt2 );

	    detail::it_algebra::substract_and_assign( x.begin() , x.end() , xtemp.begin() , xerr.begin() );
	}

	template< class DynamicalSystem , class TimeType >
        void next_step( DynamicalSystem system ,
                        ContainerType &x ,
			ContainerType &dxdt ,
                        TimeType t ,
                        TimeType dt ,
			ContainerType &xerr )
        {
	    resizer.adjust_size( x , xerr );

	    xtemp = x;
	    TimeType dt2 = 0.5*dt;

	    next_step( system , x , dxdt , t , dt );
	    next_step( system , xtemp , dxdt , t , dt2 );
	    next_step( system , xtemp , t+dt2 , dt2 );

	    detail::it_algebra::substract_and_assign( x.begin() , x.end() , xtemp.begin() , xerr.begin() );
	}


    };

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_EULER_HPP
