/*
 boost header: numeric/odeint/integrator_constant_step.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATOR_CONSTANT_STEP_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_INTEGRATOR_CONSTANT_STEP_HPP_INCLUDED

namespace boost {
namespace numeric {
namespace odeint {

    template<
	class Stepper ,
	class DynamicalSystem ,
	class TimeType ,
	class InsertIterator 
	>
    void integrate(
	Stepper stepper ,
	DynamicalSystem system ,
	TimeType dt ,
	TimeType start_time ,
	Stepper::container_type state ,
	TimeType end_time ,
	InsertIterator inserter
	)
    {
	while( start_time < end_time )
	{
	    *inserter++ = state;
	    stepper.next_step( system , state , t , dt );
	    start_time += dt;
	}
    }
    

} // odeint
} // numeric
} // boost

#endif //BOOST_NUMERIC_ODEINT_INTEGRATOR_CONSTANT_STEP_HPP_INCLUDED
