/*
 [auto_generated]
 boost/numeric/odeint/stepper/system_categories.hpp

 [begin_description]
 Definition of all system categories.
 [end_description]

 Copyright 2017 Markus Friedrich

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_SYSTEM_CATEGORIES_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_SYSTEM_CATEGORIES_HPP_INCLUDED

namespace boost {
namespace numeric {
namespace odeint {


/*
 * Tags to specify system types
 *
 * These tags can be used to detect which system category is used
 */

struct explicit_system_tag {};
struct second_order_system_tag {};
struct symplectic_system_tag {};
struct simple_symplectic_system_tag {};
struct symplectic_or_simple_symplectic_system_tag {};
struct implicit_system_tag {};


} // odeint
} // numeric
} // boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_SYSTEM_CATEGORIES_HPP_INCLUDED
