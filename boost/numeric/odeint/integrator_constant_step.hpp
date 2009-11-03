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

#include <stdexcept>

namespace boost {
namespace numeric {
namespace odeint {

    template<
	class Stepper ,
	class DynamicalSystem ,
	class TimeType ,
	class ContainerType
	>
    void integrate(
	Stepper stepper ,
	DynamicalSystem system ,
	TimeType start_time ,
	TimeType dt ,
	ContainerType &state ,
	TimeType end_time
	)
    {
	if( start_time > end_time )
	    throw std::invalid_argument( "integrate() : start_time > end_time" );

	while( start_time < end_time )
	{
	    stepper.next_step( system , state , start_time , dt );
	    start_time += dt;
	}
    }
    

} // odeint
} // numeric
} // boost

#endif //BOOST_NUMERIC_ODEINT_INTEGRATOR_CONSTANT_STEP_HPP_INCLUDED
