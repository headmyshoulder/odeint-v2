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
        class Container ,
        class Time = double ,
        class Resizer = resizer< Container >
        >
    class ode_step_euler
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
        container_type m_xtemp;
        resizer_type m_resizer;




        // public interface
    public:

        const unsigned int order() { return 1; }



        template< class DynamicalSystem >
        void next_step( DynamicalSystem system ,
                        container_type &x ,
                        const container_type &dxdt ,
                        time_type t ,
                        time_type dt )
        {
            detail::it_algebra::scale_and_add( x.begin() , x.end() , dxdt.begin() , dt );
        }



        template< class DynamicalSystem >
        void next_step( DynamicalSystem system ,
                        container_type &x ,
                        time_type t ,
                        time_type dt )
        {
            m_resizer.adjust_size( x , m_dxdt );
            system( x , m_dxdt , t );
            next_step( system , x , m_dxdt , t , dt );
        }



        template< class DynamicalSystem >
        void next_step( DynamicalSystem system ,
                        container_type &x ,
                        const container_type &dxdt ,
                        time_type t ,
                        time_type dt ,
                        container_type &xerr )
        {
            m_resizer.adjust_size( x , xerr );

            m_xtemp = x;
            time_type dt2 = 0.5 * dt;

            next_step( system , x , dxdt , t , dt );
            next_step( system , m_xtemp , dxdt , t , dt2 );
            next_step( system , m_xtemp , t+dt2 , dt2 );

            detail::it_algebra::substract_and_assign(x.begin() , 
                                                     x.end() , 
                                                     m_xtemp.begin() ,
                                                     xerr.begin() );
        }



        template< class DynamicalSystem >
        void next_step( DynamicalSystem system ,
                        container_type &x ,
                        time_type t ,
                        time_type dt ,
                        container_type &xerr )
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
