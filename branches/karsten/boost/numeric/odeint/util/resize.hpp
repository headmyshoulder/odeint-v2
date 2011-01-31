/*
 * resize.hpp
 *
 *  Created on: Jan 30, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_UTIL_RESIZE_HPP_
#define BOOST_NUMERIC_ODEINT_UTIL_RESIZE_HPP_

#include <boost/range/size.hpp>

namespace boost {
namespace numeric {
namespace odeint {


/*
 * Default implementation of resize functionality.
 * This struct has to be specialized in order to work with ublas::matrix, etc.
 */
template< class Container1 , class Container2 >
struct resize_impl
{
	static void resize( const Container1 &x1 , Container2 &x2 )
	{
		x2.resize( boost::size( x1 ) );
	}
};

template< class Container1 , class Container2 >
void resize( const Container1 &x1 , Container2 &x2 )
{
	resize_impl< Container1 , Container2 >::resize( x1 , x2 );
}



} // namespace odeint
} // namespace numeric
} // namespace boost

#endif /* BOOST_NUMERIC_ODEINT_UTIL_RESIZE_HPP_ */
