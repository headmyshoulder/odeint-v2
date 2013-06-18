/*
 [auto_generated]
 boost/numeric/odeint/external/openmp/openmp_state.hpp

 [begin_description]
 Wrappers for OpenMP.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_EXTERNAL_OPENMP_OPENMP_STATE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_EXTERNAL_OPENMP_OPENMP_STATE_HPP_INCLUDED

#include <omp.h>
#include <vector>
#include <algorithm>
#include <boost/range/adaptor/sliced.hpp>
#include <boost/numeric/odeint/util/copy.hpp>
#include <boost/numeric/odeint/util/resize.hpp>


namespace boost {
namespace numeric {
namespace odeint {

/** \brief A container that is split into distinct parts, for threading.
 */
template< class InnerState >
struct openmp_state
{
    std::vector< InnerState > m_parts;
    std::vector< size_t > offset;

    openmp_state() {}

    template< class Container >
    openmp_state( const Container &data )
    : m_parts( omp_get_num_threads() ) , offset( m_parts.size() )
    {
        const size_t part = data.size() / m_parts.size();
#       pragma omp parallel for schedule(static,1)
        for(size_t i = 0 ; i < m_parts.size() ; i++) {
            const size_t start = i * part;
            offset[i] = start;
            const size_t end = (std::min)( (i + 1) * part, data.size() );
            resize( m_parts[i], boost::adaptors::slice(data, start, end) );
            boost::numeric::odeint::copy( boost::adaptors::slice(data, start, end), m_parts[i] );
        }
    }

    InnerState & operator[](size_t i)
    {
        return m_parts[i];
    }

    const InnerState & operator[](size_t i) const
    {
        return m_parts[i];
    }

    size_t size() const
    {
        return m_parts.size();
    }

};



}
}
}


#endif

