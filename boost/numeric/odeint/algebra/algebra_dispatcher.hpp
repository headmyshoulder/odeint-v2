/*
 [auto_generated]
 boost/numeric/odeint/algebra/algebra_dispatcher.hpp

 [begin_description]
 Algebra dispatcher to automatically chose suitable algebra.
 [end_description]

 Copyright 2009-2013 Karsten Ahnert
 Copyright 2009-2013 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_NUMERIC_ODEINT_ALGEBRA_ALGEBRA_DISPATCHER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_ALGEBRA_ALGEBRA_DISPATCHER_HPP_INCLUDED


#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/array_algebra.hpp>


#include <boost/array.hpp>



namespace boost {
namespace numeric {
namespace odeint {

template< class StateType , class Enabler = void >
struct algebra_dispatcher
{
    //  range_algebra is the standard algebra
    typedef range_algebra algebra_type;
};

//specialize for array
template< class T , size_t N >
struct algebra_dispatcher< boost::array< T , N > >
{
    typedef array_algebra algebra_type;
};



}
}
}

#endif
