/*
 * destruct.hpp
 *
 *  Created on: Jan 30, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_UTIL_DESTRUCT_HPP_
#define BOOST_NUMERIC_ODEINT_UTIL_DESTRUCT_HPP_

namespace boost {
namespace numeric {
namespace odeint {

/*
 * Default implementation for destruction of a container does nothing
 * gsl_vector must be destroyed explicitly
 */
template< class Container >
struct destruct_impl
{
	static void destruct( Container &x )
	{
	}
};

template< class Container >
void destruct( Container &x )
{
	destruct_impl< Container >::destruct( x );
}


} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_UTIL_DESTRUCT_HPP_ */
