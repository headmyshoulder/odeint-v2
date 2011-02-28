/*
 * copy.hpp
 *
 *  Created on: Jan 30, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_UTIL_COPY_HPP_
#define BOOST_NUMERIC_ODEINT_UTIL_COPY_HPP_

#include <boost/range/algorithm/copy.hpp>

#include <boost/utility/enable_if.hpp>

#include <boost/numeric/odeint/util/detail/is_range.hpp>

namespace boost {
namespace numeric {
namespace odeint {

namespace detail
{
	template< class Container1 , class Container2 >
	void do_copying( const Container1 &from , Container2 &to , boost::mpl::true_ )
	{
		boost::range::copy( from , boost::begin( to ) );
	}

	template< class Container1 , class Container2 >
	void do_copying( const Container1 &from , Container2 &to , boost::mpl::false_ )
	{
		to = from;
	}
}

/*
 * Default implementation of the copy operation used the assign operator
 * gsl_vector must copied differently
 */
template< class Container1, class Container2 >
struct copy_impl
{
	static void copy( const Container1 &from , Container2 &to )
	{
		typedef typename boost::numeric::odeint::detail::is_range< Container1 >::type is_range_type;
		detail::do_copying( from , to , is_range_type() );
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
