/*
 boost header: BOOST_NUMERIC_ODEINT/thrust_resize.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_THRUST_RESIZE_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_THRUST_RESIZE_HPP_INCLUDED

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>


namespace boost {
namespace numeric {
namespace odeint {

template< class T >
struct is_resizeable< thrust::device_vector< T > >
{
	struct type : public boost::true_type { };
	const static bool value = type::value;
};

template< class T >
struct is_resizeable< thrust::host_vector< T > >
{
	struct type : public boost::true_type { };
	const static bool value = type::value;
};

} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_TR1_ARRAY_RESIZE_HPP_INCLUDED
