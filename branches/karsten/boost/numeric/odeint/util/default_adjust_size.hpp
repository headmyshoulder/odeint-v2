/*
 * adjust_size.hpp
 *
 *  Created on: Jan 30, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_UTIL_ADJUST_SIZE_HPP_
#define BOOST_NUMERIC_ODEINT_UTIL_ADJUST_SIZE_HPP_

#include <boost/numeric/odeint/util/same_size.hpp>
#include <boost/numeric/odeint/util/resize.hpp>

namespace boost {
namespace numeric {
namespace odeint {



/*
 *
 * We need this interface in order to call matrix by vector resizing,
 * i.e. mat.resize( vec.size() , vec.size() )
 *
 * See implicit_euler.hpp for usage
 */
struct default_adjust_size
{
	template< class Container1 , class Container2 >
	static bool adjust_size( const Container1 &x1 , Container2 &x2 )
	{
		if( !same_size( x1 , x2 ) )
		{
		    resize( x1 , x2 );
	        return true;
		}
		else
		{
		    return false;
		}
	}
};



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_UTIL_ADJUST_SIZE_HPP_ */
