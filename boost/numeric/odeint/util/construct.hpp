/*
 * construct.hpp
 *
 *  Created on: Jan 30, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_UTIL_CONSTRUCT_HPP_
#define BOOST_NUMERIC_ODEINT_UTIL_CONSTRUCT_HPP_

namespace boost {
namespace numeric {
namespace odeint {

/*
 * Default implementation for constructing a container does nothing
 * gsl_vector must be construct explicitly
 */
template< class Container >
struct construct_impl
{
	static void construct( Container &x )
	{
	}
};


template< class Container >
void construct( Container &x )
{
	construct_impl< Container >::construct( x );
}


} // namespace odeint
} // namespace numeric
} // namespace boost

#endif /* BOOST_NUMERIC_ODEINT_UTIL_CONSTRUCT_HPP_ */
