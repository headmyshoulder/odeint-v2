/* Boost odeint/stepper_rk4.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner
 
 This file includes the explicit 4th order runge kutta 
 solver for ordinary differential equations.

 It solves any ODE dx/dt = f(x,t).

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_RK4_HPP
#define BOOST_NUMERIC_ODEINT_STEPPER_RK4_HPP

#include <boost/concept_check.hpp>

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
    class stepper_rk4
    {

        // provide basic typedefs
    public:

        typedef Container container_type;
        typedef Resizer resizer_type;
        typedef Time time_type;
        typedef const unsigned short order_type;
        typedef typename container_type::value_type value_type;
        typedef typename container_type::iterator iterator;





        // check the concept of the ContainerType
    private:

        BOOST_CLASS_REQUIRE( container_type ,
                             boost::numeric::odeint, Container );




        // private members
    private:

        container_type m_dxdt;
        container_type m_dxt;
        container_type m_dxm;
        container_type m_dxh;
        container_type m_xt;
        resizer_type m_resizer;


	// private member functions
    private:


        // public interface
    public:

        stepper_rk4( void )
        {
        }

        order_type order() const { return 4; }

        template< class DynamicalSystem >
        void next_step( DynamicalSystem &system ,
                        container_type &x ,
                        container_type &dxdt ,
                        time_type t ,
                        time_type dt )
        {
            using namespace detail::it_algebra;

            const time_type val1 = static_cast<time_type>( 1.0 );

            m_resizer.adjust_size( x , m_dxt );
            m_resizer.adjust_size( x , m_dxm );
            m_resizer.adjust_size( x , m_xt );
            m_resizer.adjust_size( x , m_dxh );

            time_type  dh = static_cast<time_type>( 0.5 ) * dt;
            time_type th = t + dh;

            // dt * dxdt = k1
            // m_xt = x + dh*dxdt
            scale_sum( m_xt.begin(), m_xt.end(),
                       val1, x.begin(),
                       dh, dxdt.begin() );

	    // dt * m_dxt = k2
            system( m_xt , m_dxt , th );
            //m_xt = x + dh*m_dxt
            scale_sum( m_xt.begin(), m_xt.end(),
                       val1, x.begin(),
                       dh, m_dxt.begin() );

            // dt * m_dxm = k3
            system( m_xt , m_dxm , th ); 
            //m_xt = x + dt*m_dxm
            scale_sum( m_xt.begin(), m_xt.end(),
                       val1, x.begin(),
                       dt, m_dxm.begin() );

	    // dt * m_dxh = k4
            system( m_xt , m_dxh , value_type( t + dt ) );  
            //x += dt/6 * ( m_dxdt + m_dxt + val2*m_dxm )
            scale_sum( x.begin(), x.end(),
                       val1, x.begin(),
                       dt / static_cast<time_type>( 6.0 ), dxdt.begin(),
                       dt / static_cast<time_type>( 3.0 ), m_dxt.begin(),
                       dt / static_cast<time_type>( 3.0 ), m_dxm.begin(),
                       dt / static_cast<time_type>( 6.0 ), m_dxh.begin() );
        }



        template< class DynamicalSystem >
        void next_step( DynamicalSystem &system ,
                        container_type &x ,
                        time_type t ,
                        time_type dt )
        {
            m_resizer.adjust_size( x , m_dxdt );
            system( x , m_dxdt , t );
            next_step( system , x , m_dxdt , t , dt );
        }


        /* RK4 step with error estimation to 5th order according to Cash Karp */

    };

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_RK4_HPP
