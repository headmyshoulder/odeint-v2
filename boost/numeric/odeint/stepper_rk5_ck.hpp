/* Boost odeint/stepper_rk5_ck.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner
 
 This file includes the explicit 5th order runge kutta solver 
 with cash-karp error estimation (4th order) for ordinary 
 differential equations.

 It solves any ODE dx/dt = f(x,t).

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_RK5_CK_HPP
#define BOOST_NUMERIC_ODEINT_STEPPER_RK5_CK_HPP

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
    class stepper_rk5_ck
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
        container_type m_x1, m_x2, m_x3, m_x4, m_x5, m_x6;
        resizer_type m_resizer;



    public:

        order_type order() const { return 5; }

        template< class DynamicalSystem >
        void next_step( DynamicalSystem &system ,
                        container_type &x ,
                        container_type &dxdt ,
                        time_type t ,
                        time_type dt ,
                        container_type &xerr )
        {

            const time_type a2 = 0.2, a3=0.3, a4 = 0.6, a5 = 1.0, 
                a6 = 0.875, b21 = 0.2, b31 = 3.0/40.0, b32 = 9.0/40.0,
                b41 = 0.3, b42 = -0.9, b43 = 1.2, b51 = -11.0/54.0,
                b52 = 2.5, b53 = -70.0/27.0, b54 = 35.0/27.0,
                b61 = 1631.0/55296.0, b62 = 175.0/512.0, b63 = 575.0/13824.0,
                b64 = 44275.0/110592.0, b65 = 253.0/4096.0,
                c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0,
                c6 = 512.0/1771.0, dc5 = -277.0/14336.0, 
                dc1 = c1 - 2825.0/27648, dc3 = c3 - 18575.0/48384.0,
                dc4 = c4 - 13525.0/55296.0, dc6 = c6 - 0.25;
            
            using namespace detail::it_algebra;

            m_resizer.adjust_size( x , m_x1 );
            m_resizer.adjust_size( x , m_x2 );
            m_resizer.adjust_size( x , m_x3 );
            m_resizer.adjust_size( x , m_x4 );
            m_resizer.adjust_size( x , m_x5 );
            m_resizer.adjust_size( x , m_x6 );

            const time_type t_1 = static_cast<time_type>(1.0);

            //m_x1 = x + dt*b21*dxdt
            scale_sum( m_x1.begin() , m_x1.end() , 
                       t_1, x.begin() , 
                       dt*b21 , dxdt.begin() );

            system( m_x1 , m_x2 , t + dt*a2 );
            // m_x1 = x + dt*b31*dxdt + dt*b32*m_x2
            scale_sum( m_x1.begin(), m_x1.end(), 
                       t_1, x.begin(), 
                       dt*b31, dxdt.begin(),
                       dt*b32, m_x2.begin() );
            
            system( m_x1 , m_x3 , t + dt*a3 );
            // m_x1 = x + dt * (b41*dxdt + b42*m_x2 + b43*m_x3)
            scale_sum( m_x1.begin(), m_x1.end(), 
                       t_1, x.begin(),
                       dt*b41, dxdt.begin(),
                       dt*b42, m_x2.begin(),
                       dt*b43, m_x3.begin() );

            system( m_x1, m_x4 , t + dt*a4 );
            // 
            scale_sum( m_x1.begin(), m_x1.end(), 
                       t_1, x.begin(),
                       dt*b51, dxdt.begin(),
                       dt*b52, m_x2.begin(),
                       dt*b53, m_x3.begin(),
                       dt*b54, m_x4.begin() );

            system( m_x1 , m_x5 , t + dt*a5 ); // m_x5 = nr_ak5
            scale_sum( m_x1.begin(), m_x1.end(), 
                       t_1, x.begin(),
                       dt*b61, dxdt.begin(),
                       dt*b62, m_x2.begin(),
                       dt*b63, m_x3.begin(),
                       dt*b64, m_x4.begin(),
                       dt*b65, m_x5.begin() );

            system( m_x1 , m_x6 , t + dt*a6 ); // m_x6 = nr_ak6
            scale_sum( x.begin(), x.end(), 
                       t_1, x.begin(),
                       dt*c1, dxdt.begin(),
                       dt*c3, m_x3.begin(),
                       dt*c4, m_x4.begin(),
                       dt*c6, m_x6.begin() );
            
            // error estimate
            scale_sum(xerr.begin(), xerr.end(),
                      dt*dc1, dxdt.begin(),
                      dt*dc3, m_x3.begin(),
                      dt*dc4, m_x4.begin(),
                      dt*dc5, m_x5.begin(),
                      dt*dc6, m_x6.begin() );
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


}
}
}

#endif // BOOST_NUMERIC_ODEINT_STEPPER_RK5_CK_HPP
