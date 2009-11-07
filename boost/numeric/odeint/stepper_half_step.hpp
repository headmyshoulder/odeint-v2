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
	class Stepper
	>
    class ode_step_half_step
    {
        // provide basic typedefs
    public:

        typedef Stepper stepper_type;
        typedef typename Stepper::container_type container_type;
        typedef typename Stepper::resizer_type resizer_type;
        typedef typename Stepper::time_type time_type;
        typedef typename container_type::value_type value_type;
        typedef typename container_type::iterator iterator;







        // private members
    private:

        container_type m_dxdt;
        container_type m_xtemp;
        resizer_type m_resizer;
        stepper_type m_stepper;
	

        const unsigned int order() { return m_stepper.order(); }



        template< class DynamicalSystem >
        void next_step( DynamicalSystem system ,
                        container_type &x ,
                        const container_type &dxdt ,
                        time_type t ,
                        time_type dt )
        {
            m_stepper.next_step( system , x , dxdt , t , dt );
        }



        template< class DynamicalSystem >
        void next_step( DynamicalSystem system ,
                        container_type &x ,
                        time_type t ,
                        time_type dt )
        {
            m_stepper.next_step( system , x , t , dt );
        }

	/*

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

            detail::it_algebra::assign_diff(
		xerr.begin() ,
		xerr.end() ,
		x.begin() ,
		m_xtemp.begin() );
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
	*/
    }



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_EULER_HPP
