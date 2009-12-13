/* Boost odeint/stepper_midpoint.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file includes the explicit midpoint solver for
 ordinary differential equations.

 It solves any ODE dx/dt = f(x,t) via
 x(t+dt) = x(t) + dt*f(x,t)

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_MIDPOINT_HPP
#define BOOST_NUMERIC_ODEINT_STEPPER_MIDPOINT_HPP


#include <boost/numeric/odeint/detail/iterator_algebra.hpp>
#include <boost/numeric/odeint/container_traits.hpp>

#include <iostream>

namespace boost {
namespace numeric {
namespace odeint {

    template<
        class Container ,
        class Time = double ,
        class Traits = container_traits< Container >
        >
    class stepper_midpoint
    {
        // provide basic typedefs
    public:

        typedef const unsigned short order_type;
        typedef Container container_type;
        typedef Time time_type;
        typedef Traits traits_type;
        typedef typename traits_type::value_type value_type;
        typedef typename traits_type::iterator iterator;
        typedef typename traits_type::const_iterator const_iterator;




        
        // private memebers
    private:

        unsigned short m_stepcount;

        container_type m_x0;
        container_type m_x1;
        container_type m_dxdt;




    public:

        stepper_midpoint( unsigned short stepcount = 2 ) { }
        
        order_type order() const { return 2; }

        void set_stepcount( unsigned short stepcount )
        {
            if( stepcount > 1 )
                m_stepcount = stepcount;
        }

        unsigned short get_stepcount() const { return m_stepcount; }




        template< class DynamicalSystem >
        void do_step( 
                DynamicalSystem &system ,
                container_type &x ,
                container_type &dxdt ,
                time_type t ,
                time_type dt ,
                container_type &x_out )
        {
            const time_type t_1 = static_cast<time_type>( 1.0 );
            const time_type t_05 = static_cast<time_type>( 0.5 );

            const time_type h = dt/static_cast<time_type>( m_stepcount );
            const time_type h2 = static_cast<time_type>( 2.0 )*h;
            time_type th = t + h;

            traits_type::adjust_size(x, m_x0);
            traits_type::adjust_size(x, m_x1);
            traits_type::adjust_size(x, m_dxdt);

            using namespace detail::it_algebra;

            // m_x1 = x + h*dxdt
            scale_sum( traits_type::begin(m_x1), traits_type::end(m_x1),
                       t_1, traits_type::begin(x),
                       h, traits_type::begin(dxdt) );
            system( m_x1, m_dxdt, th );

            m_x0 = x;
            
            //std::clog<< "mp: " << 0 << '\t' << m_x1[0] << '\t' << m_x1[1] << '\t' << m_dxdt[0] << '\t' << m_dxdt[1] << std::endl;

            unsigned short i = 1;
            while( i != m_stepcount )
            {   // general step
                //tmp = m_x1; m_x1 = m_x0 + h2*m_dxdt; m_x0 = tmp
                scale_sum_swap( traits_type::begin(m_x1), traits_type::end(m_x1), 
                                traits_type::begin(m_x0),
                                h2, traits_type::begin(m_dxdt) );
                th += h;
                system( m_x1, m_dxdt, th);
                i++;
            }

            //std::clog<< "mp: " << i << '\t' << m_x1[0] << '\t' << m_x1[1] << '\t' << m_dxdt[0] << '\t' << m_dxdt[1] << std::endl;
            
            // last step
            // x = 0.5*( m_x0 + m_x1 + h*m_dxdt )
            scale_sum( traits_type::begin(x_out), traits_type::end(x_out),
                       t_05, traits_type::begin(m_x0),
                       t_05, traits_type::begin(m_x1),
                       t_05*h, traits_type::begin(m_dxdt) );
        }





        template< class DynamicalSystem >
        void do_step( 
                DynamicalSystem &system ,
                container_type &x ,
                container_type &dxdt ,
                time_type t ,
                time_type dt )
        {
            do_step( system, x, dxdt, t, dt, x );
        }





        template< class DynamicalSystem >
        void do_step( 
                DynamicalSystem &system ,
                container_type &x ,
                time_type t ,
                time_type dt )
        {
            traits_type::adjust_size(x, m_dxdt);
            system( x, m_dxdt, t );
            do_step( system , x, m_dxdt, t, dt );
        }
            

    };

}
}
}

#endif
