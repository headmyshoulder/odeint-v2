/*
 * same_size.hpp
 *
 *  Created on: Jan 30, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_UTIL_SAME_SIZE_HPP_
#define BOOST_NUMERIC_ODEINT_UTIL_SAME_SIZE_HPP_

#include <boost/range/size.hpp>

namespace boost {
namespace numeric {
namespace odeint {


/*
 * Default implementation of same_size functionality.
 * This struct has to be specialized in order to work with ublas::matrix, etc.
 */
template< class Container1 , class Container2 >
struct same_size_impl
{
	static bool same_size( const Container1 &x1 , const Container2 &x2 )
	{
		return ( boost::size( x1 ) == boost::size( x2 ) );
	}
};

template< class Container1 , class Container2 >
bool same_size( const Container1 &x1 , const Container2 &x2 )
{
	return same_size_impl< Container1 , Container2 >::same_size( x1 , x2 );
}



} // namespace odeint
} // namespace numeric
} // namespace boost

#endif /* BOOST_NUMERIC_ODEINT_UTIL_SAME_SIZE_HPP_ */
