/*
 * ublas_permutation_matrix_resize.hpp
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_UTIL_UBLAS_PERMUTATION_MATRIX_RESIZE_HPP_
#define BOOST_NUMERIC_ODEINT_UTIL_UBLAS_PERMUTATION_MATRIX_RESIZE_HPP_

#include <boost/type_traits/integral_constant.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace boost {
namespace numeric {
namespace odeint {


template< class T , class A >
struct is_resizeable< boost::numeric::ublas::permutation_matrix< T , A > >
{
	struct type : public boost::true_type { };
	const static bool value = type::value;
};



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_UTIL_UBLAS_PERMUTATION_MATRIX_RESIZE_HPP_ */
