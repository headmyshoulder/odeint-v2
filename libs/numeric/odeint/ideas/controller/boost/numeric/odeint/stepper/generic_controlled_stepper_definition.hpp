/*
 [auto_generated]
 boost/numeric/odeint/stepper/generic_controlled_stepper_definition.hpp

 [begin_description]
 Definition of the generic controlled stepper consisting of an error stepper, an error checker,
 and a controller.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_GENERIC_CONTROLLED_STEPPER_DEFINITION_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_GENERIC_CONTROLLED_STEPPER_DEFINITION_HPP_INCLUDED


namespace boost {
namespace numeric {
namespace odeint {


template<
    class ErrorStepper ,
    class ErrorChecker ,
    class Controller ,
    class Resizer ,
    class ErrorStepperCategory = typename ErrorStepper::stepper_category >
class generic_controlled_stepper ;


}
}
}



#endif // BOOST_NUMERIC_ODEINT_STEPPER_GENERIC_CONTROLLED_STEPPER_DEFINITION_HPP_INCLUDED
