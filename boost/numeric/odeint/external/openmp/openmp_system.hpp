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


template< class InnerState, class InnerDeriv, class Time >
struct openmp_wrapper_impl
{
    typedef boost::function<void(const InnerState &, InnerDeriv &, Time, size_t)> sys_fun_t;
    typedef openmp_state< InnerState > State;
    typedef openmp_state< InnerDeriv > Deriv;

    const sys_fun_t &f;
    openmp_wrapper_impl(const sys_fun_t &f) : f(f) {}

    inline void operator()(const State &s, Deriv &d, const Time &t, size_t off = 0 ) const
    {
#       pragma omp parallel for schedule(static,1)
        for(size_t i = 0 ; i < s.size() ; i++)
            f(s[i], d[i], t, off + s.offset[i]);
    }
};


template< class InnerState, class InnerDeriv, class Time >
inline openmp_wrapper_impl< InnerState, InnerDeriv, Time >
openmp_wrapper( const boost::function<void(const InnerState &, InnerDeriv &, Time, size_t)> &f )
{
    return openmp_wrapper_impl<InnerState, InnerDeriv, Time>(f);
}


}
}
}


#endif
