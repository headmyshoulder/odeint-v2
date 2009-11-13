/* Boost odeint/stepper_rk4.hpp header file
 
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

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_RK4_HPP
#define BOOST_NUMERIC_ODEINT_STEPPER_RK_4_HPP

#include <boost/concept_check.hpp>

#include <boost/numeric/odeint/concepts/state_concept.hpp>
#include <boost/numeric/odeint/resizer.hpp>
#include <iostream>

namespace boost {
namespace numeric {
namespace odeint {

namespace detail {
namespace constants {
    // constants for rk4 cash karp
    static const double rk4_a2 = 0.2, rk4_a3=0.3, rk4_a4 = 0.6, rk4_a5 = 1.0, 
        rk4_a6 = 0.875, rk4_b21 = 0.2, rk4_b31 = 3.0/40.0, rk4_b32 = 9.0/40.0,
        rk4_b41 = 0.3, rk4_b42 = -0.9, rk4_b43 = 1.2, rk4_b51 = -11.0/54.0,
        rk4_b52 = 2.5, rk4_b53 = -70.0/27.0, rk4_b54 = 35.0/27.0,
        rk4_b61 = 1631.0/55296.0, rk4_b62 = 175.0/512.0, rk4_b63 = 575.0/13824.0,
        rk4_b64 = 44275.0/110592.0, rk4_b65 = 253.0/4096.0,
        rk4_c1 = 37.0/378.0, rk4_c3 = 250.0/621.0, rk4_c4 = 125.0/594.0,
        rk4_c6 = 512.0/1771.0, rk4_dc5 = -277.0/14336.0, 
        rk4_dc1 = rk4_c1 - 2825.0/27648, rk4_dc3 = rk4_c3 - 18575.0/48384.0,
        rk4_dc4 = rk4_c4 - 13525.0/55296.0, rk4_dc6 = rk4_c6 - 0.25;
}
}


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
        container_type m_xt;
        container_type m_x4, m_x5, m_x6;
        resizer_type m_resizer;

        



        // public interface
    public:

        order_type order() const { return 4; }

        template< class DynamicalSystem >
        void next_step( DynamicalSystem &system ,
                        container_type &x ,
                        const container_type &dxdt ,
                        time_type t ,
                        time_type dt )
        {
            using namespace detail::it_algebra;

            const time_type val2 = time_type( 2.0 );

            m_resizer.adjust_size( x , m_dxt );
            m_resizer.adjust_size( x , m_dxm );
            m_resizer.adjust_size( x , m_xt );

            time_type  dh = time_type( 0.5 ) * dt;
            time_type th = t + dh;

            //m_xt = x + dh*dxdt
            assign_sum( m_xt.begin() , m_xt.end() , x.begin() , dxdt.begin() , dh );

            system( m_xt , m_dxt , th );
            //m_xt = x + dh*m_dxdt
            assign_sum( m_xt.begin() , m_xt.end() , x.begin() , m_dxt.begin() , dh );

            system( m_xt , m_dxm , th );
            //m_xt = x + dt*m_dxm ; m_dxm += m_dxt
            assign_sum_increment( m_xt.begin() , m_xt.end() , x.begin() ,
                                  m_dxm.begin() , m_dxt.begin() , dt );

            system( m_xt , m_dxt , value_type( t + dt ) );
            //x = dt/6 * ( m_dxdt + m_dxt + val2*m_dxm )
            increment_sum_sum( x.begin() , x.end() , dxdt.begin() ,
                               m_dxt.begin() , m_dxm.begin() ,
                               dt /  time_type( 6.0 ) , val2 );
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
        template< class DynamicalSystem >
        void next_step( DynamicalSystem &system ,
                        container_type &x ,
                        const container_type &dxdt ,
                        time_type t ,
                        time_type dt ,
                        container_type &xerr )
        {
            using namespace detail::it_algebra;
            using namespace detail::constants;
            m_resizer.adjust_size( x , m_dxt );
            m_resizer.adjust_size( x , m_dxm );
            m_resizer.adjust_size( x , m_xt );
            m_resizer.adjust_size( x , m_x4 );
            m_resizer.adjust_size( x , m_x5 );
            m_resizer.adjust_size( x , m_x6 );

            //m_xt = x + dt*b21*dxdt
            assign_sum( m_xt.begin() , m_xt.end() , x.begin() , dxdt.begin() , 
                        dt*time_type(rk4_b21) );

            system( m_xt , m_dxt , t + dt*time_type(rk4_a2) ); // m_dxt = nr_ak2
            iterator x_i = x.begin();
            iterator m_xt_i = m_xt.begin();
            typename container_type::const_iterator dxdt_i = dxdt.begin();
            iterator m_dxt_i = m_dxt.begin();
            while( x_i != x.end() ) {
                *m_xt_i++ = (*x_i++) + dt*( time_type(rk4_b31)*(*dxdt_i++) + 
                                            time_type(rk4_b32)*(*m_dxt_i++) );
            }

            system( m_xt , m_dxm , t + dt*time_type(rk4_a3) ); // m_dxm = nr_ak3
            x_i = x.begin();
            m_xt_i = m_xt.begin();
            dxdt_i = dxdt.begin();
            m_dxt_i = m_dxt.begin();
            iterator m_dxm_i = m_dxm.begin();
            while( x_i != x.end() ) {
                *m_xt_i++ = (*x_i++) + dt*( time_type(rk4_b41)*(*dxdt_i++) + 
                                            time_type(rk4_b42)*(*m_dxt_i++) + 
                                            time_type(rk4_b43)*(*m_dxm_i++) );
            }

            system( m_xt, m_x4 , t + dt*time_type(rk4_a4) ); // m_x4 = nr_ak4
            x_i = x.begin();
            m_xt_i = m_xt.begin();
            dxdt_i = dxdt.begin();
            m_dxt_i = m_dxt.begin();
            m_dxm_i = m_dxm.begin();
            iterator m_x4_i = m_x4.begin();
            while( x_i != x.end() ) {
                *m_xt_i++ = (*x_i++) + dt*( time_type(rk4_b51)*(*dxdt_i++) + 
                                            time_type(rk4_b52)*(*m_dxt_i++) +
                                            time_type(rk4_b53)*(*m_dxm_i++) + 
                                            time_type(rk4_b54)*(*m_x4_i++) ) ;
            }

            system( m_xt , m_x5 , t + dt*time_type(rk4_a5) ); // m_x5 = nr_ak5
            x_i = x.begin();
            m_xt_i = m_xt.begin();
            dxdt_i = dxdt.begin();
            m_dxt_i = m_dxt.begin();
            m_dxm_i = m_dxm.begin();
            m_x4_i = m_x4.begin();
            iterator m_x5_i = m_x5.begin();
            while( x_i != x.end() ) {
                *m_xt_i++ = (*x_i++) + dt*( time_type(rk4_b61)*(*dxdt_i++) +
                                          time_type(rk4_b62)*(*m_dxt_i++) +
                                          time_type(rk4_b63)*(*m_dxm_i++) +
                                          time_type(rk4_b64)*(*m_x4_i++) +
                                          time_type(rk4_b65)*(*m_x5_i++) );
            }

            system( m_xt , m_x6 , t + dt*time_type(rk4_a6) ); // m_x6 = nr_ak6
            x_i = x.begin();
            dxdt_i = dxdt.begin();
            m_dxm_i = m_dxm.begin();
            m_x4_i = m_x4.begin();
            iterator m_x6_i = m_x6.begin();
            while( x_i != x.end() ) {
                (*x_i++) += dt*( time_type(rk4_c1)*(*dxdt_i++) +
                                 time_type(rk4_c3)*(*m_dxm_i++) +
                                 time_type(rk4_c4)*(*m_x4_i++) +
                                 time_type(rk4_c6)*(*m_x6_i++) );
            }

            // error estimate
            iterator xerr_i = xerr.begin();
            dxdt_i = dxdt.begin();
            m_dxm_i = m_dxm.begin();
            m_x4_i = m_x4.begin();
            m_x5_i = m_x5.begin();
            m_x6_i = m_x6.begin();
            while( xerr_i != xerr.end() ) {
                *xerr_i++ = dt*( time_type(rk4_dc1)*(*dxdt_i++) +
                                  time_type(rk4_dc3)*(*m_dxm_i++) +
                                  time_type(rk4_dc4)*(*m_x4_i++) +
                                  time_type(rk4_dc5)*(*m_x5_i++) +
                                  time_type(rk4_dc6)*(*m_x6_i++) );
            }
        }

        template< class DynamicalSystem >
        void next_step( DynamicalSystem &system ,
                        container_type &x ,
                        time_type t ,
                        time_type dt ,
                        container_type &xerr )
        {
            m_resizer.adjust_size( x , m_dxdt );
            system( x , m_dxdt , t );
            next_step( system , x , m_dxdt , t , dt , xerr );
        }



    };

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_RK4_HPP
