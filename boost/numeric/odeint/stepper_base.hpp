/*
 boost header: numeric/odeint/stepper_base.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_BASE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_BASE_HPP_INCLUDED

#include <boost/noncopyable.hpp>

namespace boost {
namespace numeric {
namespace odeint {

    typedef unsigned char order_type;

    template<
	class Stepper ,
	class Container ,
	order_type Order ,
	class Time ,
	class Traits 
	>
    class explicit_stepper_base : private boost::noncopyable
    {
    public:

	// some typedef

	typedef Stepper stepper_type;
	typedef Time time_type;
	typedef Traits traits_type;
	typedef typename traits_type::container_type container_type;
        typedef typename traits_type::value_type value_type;
        typedef typename traits_type::iterator iterator;
        typedef typename traits_type::const_iterator const_iterator;




    public:

	// provide some functions

	explicit_stepper_base( stepper_type &stepper )
	    : m_stepper( stepper ) { }

	explicit_stepper_base( stepper_type &stepper , container_type &x )
	    : m_stepper( stepper )
	{
	    traits_type::adjust_size( x , m_dxdt );
	}

//	order_type order() const { return m_order; }
	const static order_type order = Order;

	template< class DynamicalSystem >
	void do_step( DynamicalSystem &system ,
		      container_type &x ,
		      time_type t ,
		      time_type dt )
	{
            system( x , m_dxdt , t );
            stepper.do_step_with_deriv( system , x , m_dxdt , t , dt );
	}

	void adjust_size( container_type &x )
	{
	    traits_type::adjust_size( x , m_dxdt );
	}

    private:

	container_type m_dxdt;
	stepper_type &m_stepper;
    };




    template<
	class ErrorStepper ,
	class Container ,
	order_type StepperOrder ,
	order_type ErrorOrder
	class Time ,
	class Traits 
	>
    class explicit_error_stepper_base : private boost::noncopyable
    {
    public:

	// some typedef

	typedef ErrorStepper stepper_type;
	typedef Time time_type;
	typedef Traits traits_type;
	typedef typename traits_type::container_type container_type;
        typedef typename traits_type::value_type value_type;
        typedef typename traits_type::iterator iterator;
        typedef typename traits_type::const_iterator const_iterator;




    public:

	// provide some functions

	explicit_stepper_base( stepper_type &stepper )
	    : m_stepper( stepper ) { }

	explicit_stepper_base( stepper_type &stepper , container_type &x )
	    : m_stepper( stepper )
	{
	    traits_type::adjust_size( x , m_dxdt );
	}

//	order_type order() const { return m_order; }
	const static order_type stepper_order = StepperOrder;
	const static order_type error_order = ErrorOrder;

	template< class DynamicalSystem >
	void do_step( DynamicalSystem &system ,
		      container_type &x ,
		      time_type t ,
		      time_type dt ,
		      container_type &err )
	{
            system( x , m_dxdt , t );
            stepper.do_step_with_deriv( system , x , m_dxdt , t , dt , err );
	}

	void adjust_size( container_type &x )
	{
	    traits_type::adjust_size( x , m_dxdt );
	}

    private:

	container_type m_dxdt;
	stepper_type &m_stepper;
    };

}
}
}

#endif //BOOST_NUMERIC_ODEINT_STEPPER_BASE_HPP_INCLUDED
