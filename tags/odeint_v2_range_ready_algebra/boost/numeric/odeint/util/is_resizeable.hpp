/*
 * is_resizeable.hpp
 *
 *  Created on: Jan 30, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_UTIL_IS_RESIZEABLE_HPP_
#define BOOST_NUMERIC_ODEINT_UTIL_IS_RESIZEABLE_HPP_

#include <vector>

#include <boost/type_traits/integral_constant.hpp>


namespace boost {
namespace numeric {
namespace odeint {

/*
 * by default any type is not resizable
 */
template< class Container >
struct is_resizeable
{
	struct type : public boost::false_type { };
	const static bool value = type::value;
};

/*
 * specialization for std::vector
 */
template< class V, class A >
struct is_resizeable< std::vector< V , A  > >
{
	struct type : public boost::true_type { };
	const static bool value = type::value;
};



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_UTIL_IS_RESIZEABLE_HPP_ */
