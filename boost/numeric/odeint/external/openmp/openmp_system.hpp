/*
 [auto_generated]
 boost/numeric/odeint/external/openmp/openmp_state.hpp

 [begin_description]
 Parallelizing wrapper for system functions
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_EXTERNAL_OPENMP_OPENMP_SYSTEM_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_EXTERNAL_OPENMP_OPENMP_SYSTEM_HPP_INCLUDED

#include "openmp_state.hpp"
#include <boost/function.hpp>


namespace boost {
namespace numeric {
namespace odeint {


template< class System >
struct openmp_wrapper_impl
{
    const System f;
    openmp_wrapper_impl(const System &f) : f(f) {}

    template< class State, class Deriv, class Time >
    inline void operator()(const State &s, Deriv &d, const Time &t, size_t off = 0 ) const
    {
#       pragma omp parallel for schedule(static,1)
        for(size_t i = 0 ; i < s.size() ; i++)
            f(s[i], d[i], t, off + s.offset[i]);
    }
};


template< class System >
inline openmp_wrapper_impl< System >
openmp_wrapper( System f )
{
    return openmp_wrapper_impl<System>(f);
}


}
}
}


#endif
