/*
 boost header: BOOST_NUMERIC_ODEINT/ublas_resize.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_UBLAS_RESIZE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_UBLAS_RESIZE_HPP_INCLUDED

#include <boost/numeric/ublas/vector.hpp>

#include <boost/type_traits/integral_constant.hpp> //for true_type and false_type

namespace boost {
namespace numeric {
namespace odeint {

/*
 * specialization for boost::numeric::ublas::vector
 */
template< class T >
struct is_resizeable< boost::numeric::ublas::vector< T > >
{
    struct type : public boost::true_type { };
    const static bool value = type::value;
};

} // odeint
} // numeric
} // boost

#endif /* BOOST_NUMERIC_ODEINT_UBLAS_RESIZE_HPP_INCLUDED */
