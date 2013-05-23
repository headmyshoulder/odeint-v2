/*
 [auto_generated]
 boost/numeric/odeint/algebra/norm_result_type.hpp

 [begin_description]
 Calculates the type of the norm_inf operation for container types
 [end_description]

 Copyright 2009-2013 Karsten Ahnert
 Copyright 2009-2013 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_NUMERIC_ODEINT_ALGEBRA_NORM_RESULT_TYPE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_ALGEBRA_NORM_RESULT_TYPE_HPP_INCLUDED

#include <boost/numeric/odeint/algebra/detail/extract_value_type.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template< typename S , typename Enabler = void >
struct norm_result_type {
    //typedef double type;
    typedef typename detail::extract_value_type< S >::type type;
};


//template< typename S >
//struct norm_result_type< S , typename S::value_type > {
//    // for anything that has value_type defined we use that
//    typedef typename detail::extract_value_type< typename S::value_type >::type type;
//};

} } }

#endif
