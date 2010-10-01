/*
   boost header: numeric/odeint/stepper_categories.hpp

   Copyright 2009 Karsten Ahnert
   Copyright 2009 Mario Mulansky

   Distributed under the Boost Software License, Version 1.0.
   (See accompanying file LICENSE_1_0.txt or
   copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_CATEGORIES_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_CATEGORIES_HPP_

namespace boost {
namespace numeric {
namespace odeint {

/*
 * Tags to specify stepper types
 */

struct explicit_error_stepper_tag {};
struct explicit_stepper_and_error_stepper_tag {};
struct explicit_error_stepper_fsal_tag {};
struct explicit_stepper_and_error_stepper_fsal_tag {};

} // odeint
} // numeric
} // boost


#endif /* BOOST_NUMERIC_ODEINT_STEPPER_CATEGORIES_HPP_ */
