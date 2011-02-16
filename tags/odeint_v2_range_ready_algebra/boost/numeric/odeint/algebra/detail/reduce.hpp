/*
 * boost header: BOOST_NUMERIC_ODEINT_ALGEBRA_DETAIL_FOR_EACH/reduce.hpp
 *
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_NUMERIC_ODEINT_ALGEBRA_DETAIL_REDUCE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_ALGEBRA_DETAIL_REDUCE_HPP_INCLUDED

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< class ValueType , class Iterator1 , class Reduction >
inline ValueType reduce( Iterator1 first1 , Iterator1 last1 , Reduction red, ValueType init)
{
	for( ; first1 != last1 ; )
		init = red( init , *first1++ );
	return init;
}

} // detail
} // odeint
} // numeric
} // boost

#endif /* BOOST_NUMERIC_ODEINT_ALGEBRA_DETAIL_REDUCE_HPP_INCLUDED */
