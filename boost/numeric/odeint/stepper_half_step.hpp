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
        // provide basic typedefs
    public:

        typedef Stepper stepper_type;
        typedef typename stepper_type::container_type container_type;
        typedef typename stepper_type::traits_type traits_type;
        typedef typename stepper_type::time_type time_type;
	typedef typename stepper_type::order_type order_type;
        typedef typename stepper_type::value_type value_type;
        typedef typename stepper_type::iterator iterator;
        typedef typename stepper_type::const_iterator const_iterator;




        // private members
    private:

        container_type m_dxdt;
        container_type m_xtemp;
        stepper_type m_stepper;
	

	// public interface
    public:

        order_type order() const { return m_stepper.order(); }

        order_type order_error() const 
        {   /* Order of the error term is the order of the underlying stepper + 1 */
            return m_stepper.order() + 1; 
        }

        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                        container_type &x ,
                        container_type &dxdt ,
                        time_type t ,
                        time_type dt )
        {
            m_stepper.do_step( system , x , dxdt , t , dt );
        }



        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                        container_type &x ,
                        time_type t ,
                        time_type dt )
        {
            m_stepper.do_step( system , x , t , dt );
        }

        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                        container_type &x ,
                        container_type &dxdt ,
                        time_type t ,
                        time_type dt ,
                        container_type &xerr )
        {
            traits_type::adjust_size( x , xerr );

            /*** BUG FOR BLITZ ARRAYS ***/
            /* for blitzz arrays, the copy constructor creates a REFERENCE and not a copy!!! */
            traits_type::adjust_size( x , m_xtemp );
            m_xtemp = container_type(x);
            /*** END BUG ***/

            time_type dt2 = static_cast<time_type>(0.5) * dt;

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



        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                        container_type &x ,
                        time_type t ,
                        time_type dt ,
                        container_type &xerr )
        {
            traits_type::adjust_size( x , m_dxdt );
            system( x , m_dxdt , t );
            do_step( system , x , m_dxdt , t , dt , xerr );
        }
    };



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_HALF_STEP_HPP
