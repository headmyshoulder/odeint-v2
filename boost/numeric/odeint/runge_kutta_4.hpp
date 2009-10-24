/* Boost odeint/runge_kutta_4.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner
 
 This file includes the explicit runge kutta solver for
 ordinary differential equations.

 It solves any ODE dx/dt = f(x,t).

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_RUNGE_KUTTA_4_HPP
#define BOOST_NUMERIC_ODEINT_RUNGE_KUTTA_4_HPP

#include <boost/concept_check.hpp>

#include <boost/numeric/odeint/concepts/state_concept.hpp>
#include <boost/numeric/odeint/resizer.hpp>

namespace boost {
namespace numeric {
namespace odeint {


    template<
	class ContainerType ,
	class ResizeType = resizer< ContainerType >
	>
    class ode_step_runge_kutta_4
    {
        BOOST_CLASS_REQUIRE( ContainerType , boost::numeric::odeint, StateType );

        ContainerType dxdt;
        ContainerType dxt;
        ContainerType dxm;
	ContainerType xt;

        ResizeType resizer;

        typedef typename ContainerType::iterator iterator;
	typedef typename ContainerType::value_type value_type;
        
    public:

        template< class DynamicalSystem , class TimeType>
        void next_step( DynamicalSystem system ,
                        ContainerType &x ,
                        TimeType t ,
                        TimeType dt )
        {
	    const value_type val2 = value_type( 2.0 );

            if( ! resizer.same_size( x , dxdt ) ) resizer.resize( x , dxdt );
            if( ! resizer.same_size( x , dxt ) ) resizer.resize( x , dxt );
            if( ! resizer.same_size( x , dxm ) ) resizer.resize( x , dxm );
            if( ! resizer.same_size( x , xt ) ) resizer.resize( x , xt );

	    value_type dt_val = value_type( dt );
            value_type dh = dt_val * value_type( 0.5 );
	    value_type d6 = dt_val / value_type( 6.0 );
            value_type th = dh + value_type( t );

	    iterator iter1 , iter2 ,iter3 , iter4;
	    iterator x_end = x.end() , xt_end = xt.end();

            system( x , dxdt , t );
	    iter1 = xt.begin() ; iter2 = x.begin() ; iter3 = dxdt.begin();
	    while( iter1 != xt_end )
		(*iter1++) = (*iter2++) + dh * (*iter3++);

	    system( xt , dxt , th );
	    iter1 = xt.begin() ; iter2 = x.begin() ; iter3 = dxt.begin();
	    while( iter1 != xt_end ) 
		(*iter1++) = (*iter2++) + dh * (*iter3++);

	    system( xt , dxm , th );
	    iter1 = xt.begin() ; iter2 = x.begin() ; iter3 = dxm.begin() ; iter4  = dxt.begin();
	    while( iter1 != xt_end )
	    {
		(*iter1++) = (*iter2++) + dt_val * (*iter3);
		(*iter3++) += (*iter4++);
	    }

	    system( xt , dxt , value_type( t + dt ) );
	    iter1 = x.begin() ; iter2 = dxdt.begin() ; iter3 = dxt.begin() ; iter4 = dxm.begin();
	    while( iter1 != x_end )
		(*iter1++) += d6 * ( (*iter2++) + (*iter3++) + val2 * (*iter4++) );
        }
    };

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_RUNGE_KUTTA_4_HPP
