/*
 [auto_generated]
 boost/numeric/odeint/stepper/generic_controlled_stepper.hpp

 [begin_description]
 tba.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_GENERIC_CONTROLLED_STEPPER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_GENERIC_CONTROLLED_STEPPER_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

namespace boost {
namespace numeric {
namespace odeint {


template<
    class ErrorStepper ,
    class ErrorChecker ,
    class Controller ,
    class Resizer ,
    class ErrorStepperCategorie = typename ErrorStepper::stepper_type >
class generic_controlled_stepper ;



}
}
}



#endif // BOOST_NUMERIC_ODEINT_STEPPER_GENERIC_CONTROLLED_STEPPER_HPP_INCLUDED
