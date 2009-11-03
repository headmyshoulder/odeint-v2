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

        ContainerType m_dxdt;
	ContainerType m_xtemp;
        ResizeType m_resizer;

    public:


	// provide ContainerType, ResizeType, iterator and value_type to users of this class
	typedef ContainerType container_type;
	typedef ResizerType resizer_type;
        typedef typename container_type::iterator iterator;
	typedef typename container_type::value_type value_type;

	const unsigned int order() { return 1; }

	template< class DynamicalSystem , class TimeType >
	void next_step( DynamicalSystem system ,
                        ContainerType &x ,
			const ContainerType &dxdt ,
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
	    m_resizer.adjust_size( x , m_dxdt );
            system( x , m_dxdt , t );
	    next_step( system , x , m_dxdt , t , dt );
        }

	template< class DynamicalSystem , class TimeType >
	void next_step( DynamicalSystem system ,
			ContainerType &x ,
			const ContainerType &dxdt ,
			TimeType t ,
			TimeType dt ,
			ContainerType &xerr )
        {
	    m_resizer.adjust_size( x , xerr );

	    m_xtemp = x;
	    TimeType dt2 = 0.5 * dt;

	    next_step( system , x , dxdt , t , dt );
	    next_step( system , m_xtemp , dxdt , t , dt2 );
	    next_step( system , m_xtemp , t+dt2 , dt2 );

	    detail::it_algebra::substract_and_assign( x.begin() , x.end() , m_xtemp.begin() , xerr.begin() );
	}





	template< class DynamicalSystem , class TimeType >
        void next_step( DynamicalSystem system ,
                        ContainerType &x ,
                        TimeType t ,
                        TimeType dt ,
			ContainerType &xerr )
        {
	    m_resizer.check_size_and_resize( x , m_dxdt );
	    system( x , m_dxdt , t );
	    next_step( system , x , m_dxdt , t , dt , xerr );
	}
    };

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_EULER_HPP
