/* Boost odeint/stepper_rk5_ck.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner
 
 This file includes the explicit 5th order runge kutta solver 
 with cash-karp error estimation (4th order) for ordinary 
 differential equations.

 doi:10.1145/79505.79507.
 http://portal.acm.org/citation.cfm?doid=79505.79507

 It solves any ODE dx/dt = f(x,t).

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_RK5_CK_HPP
#define BOOST_NUMERIC_ODEINT_STEPPER_RK5_CK_HPP

#include <boost/numeric/odeint/concepts/state_concept.hpp>

namespace boost {
namespace numeric {
namespace odeint {

    template<
        class Container ,
        class Time = double ,
        class Traits = container_traits< Container >
        >
    class stepper_rk5_ck
    {

        // provide basic typedefs
    public:

        typedef unsigned short order_type;
        typedef Time time_type;
        typedef Traits traits_type;
        typedef typename traits_type::container_type container_type;
        typedef typename traits_type::value_type value_type;
//        typedef typename traits_type::iterator iterator;
//        typedef typename traits_type::const_iterator const_iterator;





        // private members
    private:

        container_type m_dxdt;
        container_type m_x1, m_x2, m_x3, m_x4, m_x5, m_x6;



    public:

	// the 5th order evaluation is used as step result
        order_type order_error_step() const { return 5; }


        order_type order_error() const { return 5; }


	// standard constructor, the internal container are uninitialized
	stepper_rk5_ck( void )
	{
	}


	// contructor which adjusts the internal containers
	stepper_rk5_ck( const container_type &x )
	{
	    adjust_size( x );
	}


	void adjust_size( const container_type &x )
	{
            traits_type::adjust_size( x , m_dxdt );
            traits_type::adjust_size( x , m_x1 );
            traits_type::adjust_size( x , m_x2 );
            traits_type::adjust_size( x , m_x3 );
            traits_type::adjust_size( x , m_x4 );
            traits_type::adjust_size( x , m_x5 );
            traits_type::adjust_size( x , m_x6 );
	}




        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                      container_type &x ,
                      const container_type &dxdt ,
                      time_type t ,
                      time_type dt ,
                      container_type &xerr )
        {

            const time_type a2 = static_cast<time_type> ( 0.2 );
	    const time_type a3 = static_cast<time_type> ( 0.3 );
	    const time_type a4 = static_cast<time_type> ( 0.6 );
	    const time_type a5 = static_cast<time_type> ( 1.0 );
	    const time_type a6 = static_cast<time_type> ( 0.875 );

	    const time_type b21 = static_cast<time_type> ( 0.2 );
	    const time_type b31 = static_cast<time_type> ( 3.0 ) / static_cast<time_type>( 40.0 );
	    const time_type b32 = static_cast<time_type> ( 9.0 ) / static_cast<time_type>( 40.0 );
	    const time_type b41 = static_cast<time_type> ( 0.3 );
	    const time_type b42 = static_cast<time_type> ( -0.9 );
	    const time_type b43 = static_cast<time_type> ( 1.2 );
	    const time_type b51 = static_cast<time_type> ( -11.0 ) / static_cast<time_type>( 54.0 );
	    const time_type b52 = static_cast<time_type> ( 2.5 );
	    const time_type b53 = static_cast<time_type> ( -70.0 ) / static_cast<time_type>( 27.0 );
	    const time_type b54 = static_cast<time_type> ( 35.0 ) / static_cast<time_type>( 27.0 );
	    const time_type b61 = static_cast<time_type> ( 1631.0 ) / static_cast<time_type>( 55296.0 );
	    const time_type b62 = static_cast<time_type> ( 175.0 ) / static_cast<time_type>( 512.0 );
	    const time_type b63 = static_cast<time_type> ( 575.0 ) / static_cast<time_type>( 13824.0 );
	    const time_type b64 = static_cast<time_type> ( 44275.0 ) / static_cast<time_type>( 110592.0 );
	    const time_type b65 = static_cast<time_type> ( 253.0 ) / static_cast<time_type>( 4096.0 );

	    const time_type c1 = static_cast<time_type> ( 37.0 ) / static_cast<time_type>( 378.0 );
	    const time_type c3 = static_cast<time_type> ( 250.0 ) / static_cast<time_type>( 621.0 );
	    const time_type c4 = static_cast<time_type> ( 125.0 ) / static_cast<time_type>( 594.0 );
	    const time_type c6 = static_cast<time_type> ( 512.0 ) / static_cast<time_type>( 1771.0 );

	    const time_type dc1 = c1 - static_cast<time_type> ( 2825.0 ) / static_cast<time_type>( 27648 );
	    const time_type dc3 = c3 - static_cast<time_type> ( 18575.0 ) / static_cast<time_type>( 48384.0 );
	    const time_type dc4 = c4 - static_cast<time_type> ( 13525.0 ) / static_cast<time_type>( 55296.0 );
	    const time_type dc5 = static_cast<time_type> ( -277.0 ) / static_cast<time_type>( 14336.0 );
	    const time_type dc6 = c6 - static_cast<time_type> ( 0.25 );
            
            using namespace detail::it_algebra;


            const time_type t_1 = static_cast<time_type>(1.0);

            //m_x1 = x + dt*b21*dxdt
            scale_sum( traits_type::begin(m_x1) , traits_type::end(m_x1) , 
                       t_1, traits_type::begin(x) , 
                       dt*b21 , traits_type::begin(dxdt) );


            system( m_x1 , m_x2 , t + dt*a2 );
            // m_x1 = x + dt*b31*dxdt + dt*b32*m_x2
            scale_sum( traits_type::begin(m_x1), traits_type::end(m_x1), 
                       t_1, traits_type::begin(x), 
                       dt*b31, traits_type::begin(dxdt),
                       dt*b32, traits_type::begin(m_x2) );
            

            system( m_x1 , m_x3 , t + dt*a3 );
            // m_x1 = x + dt * (b41*dxdt + b42*m_x2 + b43*m_x3)
            scale_sum( traits_type::begin(m_x1), traits_type::end(m_x1), 
                       t_1, traits_type::begin(x),
                       dt*b41, traits_type::begin(dxdt),
                       dt*b42, traits_type::begin(m_x2),
                       dt*b43, traits_type::begin(m_x3) );


            system( m_x1, m_x4 , t + dt*a4 );
            // 
            scale_sum( traits_type::begin(m_x1), traits_type::end(m_x1), 
                       t_1, traits_type::begin(x),
                       dt*b51, traits_type::begin(dxdt),
                       dt*b52, traits_type::begin(m_x2),
                       dt*b53, traits_type::begin(m_x3),
                       dt*b54, traits_type::begin(m_x4) );


            system( m_x1 , m_x5 , t + dt*a5 ); // m_x5 = nr_ak5
            scale_sum( traits_type::begin(m_x1), traits_type::end(m_x1), 
                       t_1, traits_type::begin(x),
                       dt*b61, traits_type::begin(dxdt),
                       dt*b62, traits_type::begin(m_x2),
                       dt*b63, traits_type::begin(m_x3),
                       dt*b64, traits_type::begin(m_x4),
                       dt*b65, traits_type::begin(m_x5) );


            system( m_x1 , m_x6 , t + dt*a6 ); // m_x6 = nr_ak6
            scale_sum( traits_type::begin(x), traits_type::end(x), 
                       t_1, traits_type::begin(x),
                       dt*c1, traits_type::begin(dxdt),
                       dt*c3, traits_type::begin(m_x3),
                       dt*c4, traits_type::begin(m_x4),
                       dt*c6, traits_type::begin(m_x6) );
            

            // error estimate
            scale_sum(traits_type::begin(xerr), traits_type::end(xerr),
                      dt*dc1, traits_type::begin(dxdt),
                      dt*dc3, traits_type::begin(m_x3),
                      dt*dc4, traits_type::begin(m_x4),
                      dt*dc5, traits_type::begin(m_x5),
                      dt*dc6, traits_type::begin(m_x6) );
        }




        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                        container_type &x ,
                        time_type t ,
                        time_type dt ,
                        container_type &xerr )
        {
            system( x , m_dxdt , t );
            do_step( system , x , m_dxdt , t , dt , xerr );
        }


    };


}
}
}

#endif // BOOST_NUMERIC_ODEINT_STEPPER_RK5_CK_HPP
