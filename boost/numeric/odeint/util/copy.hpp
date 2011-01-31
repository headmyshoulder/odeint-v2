/*
 * copy.hpp
 *
 *  Created on: Jan 30, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_UTIL_COPY_HPP_
#define BOOST_NUMERIC_ODEINT_UTIL_COPY_HPP_

namespace boost {
namespace numeric {
namespace odeint {

/*
 * Default implementation of the copy operation used the assign operator
 * gsl_vector must copied differently
 */
template< class Container1, class Container2 >
struct copy_impl
{
	static void copy( const Container1 &from , Container2 &to )
	{
		to = from;
	}
};

template< class Container1 , class Container2 >
void copy( const Container1 &from , Container2 &to )
{
	copy_impl< Container1 , Container2 >::copy( from , to );
}


} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_UTIL_COPY_HPP_ */
