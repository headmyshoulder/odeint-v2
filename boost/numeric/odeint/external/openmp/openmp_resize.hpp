/*
 [auto_generated]
 boost/numeric/odeint/external/openmp/openmp_resize.hpp

 [begin_description]
 Delegate resizing of OpenMP-state to inner state.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_EXTERNAL_OPENMP_OPENMP_RESIZE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_EXTERNAL_OPENMP_OPENMP_RESIZE_HPP_INCLUDED

#include <boost/numeric/odeint/util/resizer.hpp>
#include "openmp_state.hpp"

namespace boost {
namespace numeric {
namespace odeint {


template< class InnerState >
struct is_resizeable< openmp_state< InnerState > >
{
    typedef typename is_resizeable< InnerState >::type type;
    const static bool value = is_resizeable< InnerState >::value;
};


template< class InnerState1, class InnerState2 >
struct same_size_impl< openmp_state< InnerState1 > , openmp_state< InnerState2 > >
{
    static bool same_size( const openmp_state< InnerState1 > &x , const openmp_state< InnerState2 > &y )
    {
        if( x.size() != y.size() ) return false;
        for( size_t i = 0 ; i != x.size() ; i++ )
            if( !same_size(x[i], y[i]) ) return false;
        return true;
    }
};


template< class InnerStateOut, class InnerStateIn >
struct resize_impl< openmp_state< InnerStateOut > , openmp_state< InnerStateIn > >
{
    static void resize( openmp_state< InnerStateOut > &x , const openmp_state< InnerStateIn > &y )
    {
        x.m_parts.resize( y.size() );
        x.offset.resize( y.size() );
#       pragma omp parallel for schedule(static,1)
        for(size_t i = 0 ; i < x.size() ; i++)
            boost::numeric::odeint::resize( x[i], y[i] );
        size_t off = 0;
        for(size_t i = 0 ; i < x.size() ; i++) {
            x.offset[i] = off;
            off += x[i].size();
        }
    }
};


template< class InnerStateIn , class InnerStateOut >
struct copy_impl< openmp_state< InnerStateIn >, openmp_state< InnerStateOut > >
{
    static void copy( const openmp_state< InnerStateIn > &from, openmp_state< InnerStateOut > &to )
    {
#       pragma omp parallel for schedule(static,1)
        for(size_t i = 0 ; i < from.size() ; i++)
            copy( from[i], to[i] );
    }
};


} // odeint
} // numeric
} // boost


#endif // BOOST_NUMERIC_ODEINT_EXTERNAL_THRUST_THRUST_RESIZE_HPP_INCLUDED

