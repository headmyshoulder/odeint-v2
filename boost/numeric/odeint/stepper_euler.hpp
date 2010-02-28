/* Boost odeint/stepper_euler.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner
 
 This file includes the explicit euler solver for
 ordinary differential equations.

 It solves any ODE dx/dt = f(x,t) via
 x(t+dt) = x(t) + dt*f(x,t)

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_EULER_HPP
#define BOOST_NUMERIC_ODEINT_STEPPER_EULER_HPP

#include <boost/numeric/odeint/detail/iterator_algebra.hpp>
#include <boost/numeric/odeint/container_traits.hpp>

#include <boost/numeric/odeint/stepper_base.hpp>


namespace boost {
namespace numeric {
namespace odeint {

    template<
        class Container ,
        class Time = double ,
        class Traits = container_traits< Container >
        >
    class stepper_euler : public stepper_base<
	stepper_euler< Container , Time , Traits > ,
	Container ,
	1 ,
	Time ,
	Traits >
    {
        // provide basic typedefs
    public:

	typedef stepper_base< stepper_euler< Container , Time , Traits > ,
			      Container , 1 , Time , Traits > base_type;

        typedef typename base_type::time_type time_type;
//        typedef typename base_type::order_type order_type;
/*
        typedef unsigned short order_type;
        typedef Time time_type;
*/

        typedef Traits traits_type;
        typedef typename traits_type::container_type container_type;
        typedef typename traits_type::value_type value_type;
        typedef typename traits_type::iterator iterator;
        typedef typename traits_type::const_iterator const_iterator;



        // private members
    private:

        container_type m_dxdt;



        // public interface
    public:

//        order_type order() const { return 1; }


        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
		      container_type &x ,
		      const container_type &dxdt ,
		      time_type t ,
		      time_type dt )
        {
            //x = x + dt*dxdt
            detail::it_algebra::increment( traits_type::begin(x) ,
                                           traits_type::end(x) ,
                                           traits_type::begin(dxdt) , 
                                           dt );
        }

/*	template< class DynamicalSystem >
	void do_step( DynamicalSystem &system ,
		      container_type &x ,
		      time_type t ,
		      time_type dt )
	{
            traits_type::adjust_size( x , m_dxdt );
            system( x , m_dxdt , t );
            do_step( system , x , m_dxdt , t , dt );
	    }*/
    };



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_EULER_HPP
