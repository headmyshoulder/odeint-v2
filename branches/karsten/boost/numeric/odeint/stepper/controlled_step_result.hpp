/*
 * controlled_step_result.hpp
 *
 *  Created on: Jan 27, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_STEPPER_RESULT_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_STEPPER_RESULT_HPP_


namespace boost {
namespace numeric {
namespace odeint {

typedef enum
{
    success_step_size_unchanged ,
    step_size_decreased ,
    success_step_size_increased
} controlled_step_result;

} // namespace odeint
} // numeric
} // boost

#endif /* BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_STEPPER_RESULT_HPP_ */
