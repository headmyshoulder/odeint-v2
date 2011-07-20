/*
 [auto_generated]
 boost/numeric/odeint/external/thrust/thrust_resize.hpp

 [begin_description]
 Enable resizing for thrusts device and host_vector.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_EXTERNAL_THRUST_THRUST_RESIZE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_EXTERNAL_THRUST_THRUST_RESIZE_HPP_INCLUDED


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


#endif // BOOST_NUMERIC_ODEINT_EXTERNAL_THRUST_THRUST_RESIZE_HPP_INCLUDED
