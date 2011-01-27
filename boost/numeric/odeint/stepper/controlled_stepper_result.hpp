/*
 * controlled_stepper_result.hpp
 *
 *  Created on: Jan 27, 2011
 *      Author: karsten
 */

#ifndef CONTROLLED_STEPPER_RESULT_HPP_
#define CONTROLLED_STEPPER_RESULT_HPP_


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

#endif /* CONTROLLED_STEPPER_RESULT_HPP_ */
