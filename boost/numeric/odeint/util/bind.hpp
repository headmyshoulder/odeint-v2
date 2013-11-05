/*
 *     [begin_description]
 *     Boost bind pull the placeholders, _1, _2, ... into global
 *     namespace. This can conflict with the C++03 TR1 and C++11 
 *     std::placeholders. This header provides a workaround for 
 *     this problem.
 *     [end_description]
 *        
 *     Copyright 2012 Christoph Koke
 *           
 *     Distributed under the Boost Software License, Version 1.0.
 *     (See accompanying file LICENSE_1_0.txt or
 *     copy at http://www.boost.org/LICENSE_1_0.txt)
 * */

#ifndef BOOST_NUMERIC_ODEINT_UTIL_BIND_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_UTIL_BIND_HPP_INCLUDED


#include <boost/config.hpp>


#ifdef BOOST_NO_CXX11_HDR_FUNCTIONAL
#include <boost/bind.hpp>
#else
#include <functional>
#endif


namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

    
#ifdef BOOST_NO_CXX11_HDR_FUNCTIONAL

using ::boost::bind;
using ::_1;
using ::_2;

#else

using ::std::bind;
using namespace ::std::placeholders;

#endif


}
}
}
}

#endif // BOOST_NUMERIC_ODEINT_UTIL_BIND_HPP_INCLUDED
