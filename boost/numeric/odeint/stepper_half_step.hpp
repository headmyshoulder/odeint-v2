/* Boost odeint/stepper_half_stepr.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner
 
 This file includes a stepper which calculates the
 error during one step from performing two steps with
 the halt stepsize. It works with arbitray steppers

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_HALF_STEP_HPP
#define BOOST_NUMERIC_ODEINT_STEPPER_HALF_STEP_HPP

#include <boost/numeric/odeint/detail/iterator_algebra.hpp>

#include <iostream>

namespace boost {
namespace numeric {
namespace odeint {


    template< class Stepper >
    class stepper_half_step
    {
        //
        // provide basic typedefs
        //
    public:

        typedef Stepper stepper_type;
        typedef typename stepper_type::container_type container_type;
        typedef typename stepper_type::traits_type traits_type;
        typedef typename stepper_type::time_type time_type;
	typedef typename stepper_type::order_type order_type;
        typedef typename stepper_type::value_type value_type;
        typedef typename stepper_type::iterator iterator;
        typedef typename stepper_type::const_iterator const_iterator;



        //
        // private members
        //
    private:

        container_type m_dxdt;
        container_type m_xtemp;
        stepper_type m_stepper;
	

        //
	// public interface
        //
    public:


        // standard constructor
        stepper_half_step( void )
        {
        }

        // contructor, which adjust the size of internal containers
        stepper_half_step( const container_type &x )
        {
            adjust_size( x );
        }

        // adjust the size of m_dxdt , m_xtemp und m_stepper
        void adjust_size( const container_type &x )
        {
            m_stepper.adjust_size( x );
            traits_type::adjust_size( x , m_dxdt );
            traits_type::adjust_size( x , m_xtemp );
        }

        // the order of the step if a normal step is performed
        order_type order_step( void ) const
        {

            return m_stepper.order_step();
        }

        // the order of the step if an error step is performed
        order_type order_error_step( void ) const
        {

            return m_stepper.order_step();
        }

        // the order of the error term if the error step is performed
        order_type order_error( void ) const 
        {

            return m_stepper.order_step() + 1; 
        }


        // performs a normal step, without error calculation
        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                      container_type &x ,
                      const container_type &dxdt ,
                      time_type t ,
                      time_type dt )
        {
            m_stepper.do_step( system , x , dxdt , t , dt );
        }


        // performs a normal step, without error calculation
        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                      container_type &x ,
                      time_type t ,
                      time_type dt )
        {
            m_stepper.do_step( system , x , t , dt );
        }



        // performs a error step with error calculation
        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                        container_type &x ,
                        container_type &dxdt ,
                        time_type t ,
                        time_type dt ,
                        container_type &xerr )
        {
            time_type dt2 = static_cast<time_type>(0.5) * dt;

            m_xtemp = x;

            do_step( system , m_xtemp , dxdt , t , dt );
            do_step( system , x , dxdt , t , dt2 );
            do_step( system , x , t+dt2 , dt2 );

            detail::it_algebra::scale_sum( traits_type::begin(xerr) ,
                                           traits_type::end(xerr) ,
                                           static_cast< value_type >(1.0),
                                           traits_type::begin(m_xtemp) ,
                                           static_cast< value_type >(-1.0),
                                           traits_type::begin(x) );
        }




        // performs a error step with error calculation
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



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_HALF_STEP_HPP
