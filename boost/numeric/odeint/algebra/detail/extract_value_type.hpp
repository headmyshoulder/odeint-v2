/*
 [auto_generated]
 boost/numeric/odeint/algebra/detail/extract_value_type.hpp

 [begin_description]
 Extract true value type from complex types (eg. std::complex)
 [end_description]

 Copyright 2009-2013 Karsten Ahnert
 Copyright 2009-2013 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_NUMERIC_ODEINT_ALGEBRA_DETAIL_EXTRACT_VALUE_TYPE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_ALGEBRA_DETAIL_EXTRACT_VALUE_TYPE_HPP_INCLUDED

#include <complex>

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< typename S >
struct extract_value_type { typedef S type; };

template< typename T >
struct extract_value_type< std::complex< T > > { typedef T type; };

} } } }

#endif
