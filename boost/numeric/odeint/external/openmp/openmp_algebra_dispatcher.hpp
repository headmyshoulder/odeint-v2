/*
 [auto_generated]
 boost/numeric/odeint/external/openmp/openmp_algebra_dispatcher.hpp

 [begin_description]
 Algebra dispatcher to automatically chose suitable algebra.
 [end_description]

 Copyright 2009-2013 Karsten Ahnert
 Copyright 2009-2013 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_NUMERIC_ODEINT_OPENMP_OPENMP_ALGEBRA_DISPATCHER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_OPENMP_OPENMP_ALGEBRA_DISPATCHER_HPP_INCLUDED

#include <boost/numeric/odeint/external/openmp/openmp_algebra.hpp>
#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include "openmp_state.hpp"

namespace boost {
namespace numeric {
namespace odeint {

template< class Inner >
struct algebra_dispatcher< openmp_state< Inner > >
{
    typedef openmp_algebra<
        typename algebra_dispatcher< Inner >::algebra_type
    > algebra_type;
};

}
}
}

#endif
