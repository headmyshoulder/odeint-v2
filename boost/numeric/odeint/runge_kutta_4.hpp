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
#include <iostream>

namespace boost {
namespace numeric {
namespace odeint {


    template<
        class Container ,
	class Time = double ,
        class Resizer = resizer< Container >
        >
    class ode_step_runge_kutta_4
    {

        // provide basic typedefs
    public:

        typedef Container container_type;
        typedef Resizer resizer_type;
        typedef Time time_type;
        typedef typename container_type::value_type value_type;
        typedef typename container_type::iterator iterator;





        // check the concept of the ContainerType
    private:

        BOOST_CLASS_REQUIRE( container_type ,
			     boost::numeric::odeint, StateType );




        // private members
    private:

        container_type m_dxdt;
        container_type m_dxt;
        container_type m_dxm;
        container_type m_xt;
        resizer_type m_resizer;

        



	// public interface
    public:

        template< class DynamicalSystem >
        void next_step( DynamicalSystem system ,
                        container_type &x ,
                        time_type t ,
                        time_type dt )
        {
	    using namespace detail::it_algebra;

            const time_type val2 = time_type( 2.0 );

	    m_resizer.adjust_size( x , m_dxdt );
	    m_resizer.adjust_size( x , m_dxt );
	    m_resizer.adjust_size( x , m_dxm );
	    m_resizer.adjust_size( x , m_xt );

            time_type  dh = time_type( 0.5 ) * dt;
            time_type d6 = dt /  time_type( 6.0 );
            time_type th = t + dh;

            iterator iter1 , iter2 ,iter3 , iter4;
            iterator x_end = x.end() , xt_end = m_xt.end();

            system( x , m_dxdt , t );
	    scale_and_add_and_assign( x.begin() , x.end() , m_dxdt.begin() , m_xt.begin() , dh );

            system( m_xt , m_dxt , th );
	    scale_and_add_and_assign( x.begin() , x.end() , m_dxt.begin() , m_xt.begin() , dh );

            system( m_xt , m_dxm , th );
            iter1 = m_xt.begin() ; iter2 = x.begin() ; iter3 = m_dxm.begin() ; iter4  = m_dxt.begin();
            while( iter1 != xt_end )
            {
                (*iter1++) = (*iter2++) + dt * (*iter3);
                (*iter3++) += (*iter4++);
            }

            system( m_xt , m_dxt , value_type( t + dt ) );
	    scale_and_add_and_add_and_assign( m_dxdt.begin() , m_dxdt.end() , m_dxt.begin() , m_dxm.begin() , x.begin() , val2 , d6 );

        }


    };

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_RUNGE_KUTTA_4_HPP
