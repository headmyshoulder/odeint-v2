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

#if __cplusplus >= 201103L
#define BOOST_NUMERIC_ODEINT_CXX11 1
#endif

#if BOOST_NUMERIC_ODEINT_CXX11 
#include <functional>
#include <type_traits>
#else
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#endif

namespace boost {

#if BOOST_NUMERIC_ODEINT_CXX11
template<typename T> class reference_wrapper;

template<typename T> class unwrap_reference;
#endif

namespace numeric {
namespace odeint {
namespace detail {

#if BOOST_NUMERIC_ODEINT_CXX11 
using ::std::bind;
using ::std::ref;
using namespace ::std::placeholders;


template<typename T>
struct unwrap_reference
{
	typedef typename std::remove_reference<T>::type type;
};

template<typename T>
struct unwrap_reference< std::reference_wrapper<T> > 
{
	typedef typename std::remove_reference<T>::type type;
};

template<typename T>
struct unwrap_reference< boost::reference_wrapper<T> >
{
	    typedef typename boost::unwrap_reference<T>::type type;
};

#else
using ::boost::bind;
using ::boost::ref;
using ::boost::unwrap_reference;
using ::_1;
using ::_2;
#endif

}
}
}
}

#endif // BOOST_NUMERIC_ODEINT_UTIL_BIND_HPP_INCLUDED
