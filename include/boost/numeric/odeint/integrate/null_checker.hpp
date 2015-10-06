/*
 [auto_generated]
 boost/numeric/odeint/integrate/null_checker.hpp

 [begin_description]
 null_checker
 [end_description]

 Copyright 2015 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_NULL_CHECKER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_INTEGRATE_NULL_CHECKER_HPP_INCLUDED

namespace boost {
namespace numeric {
namespace odeint {

struct null_checker
{
    void operator()( void ) const
    { }

    void reset( void ) const
    { }
};

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_ODEINT_INTEGRATE_NULL_CHECKER_HPP_INCLUDED
