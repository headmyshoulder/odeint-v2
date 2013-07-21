/*
 [auto_generated]
 boost/numeric/odeint/external/mpi/mpi_vector_state.hpp

 [begin_description]
 Copying a container from/to an mpi_state splits/joins it.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_EXTERNAL_MPI_MPI_VECTOR_STATE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_EXTERNAL_MPI_MPI_VECTOR_STATE_HPP_INCLUDED

#include <vector>
#include <algorithm>
#include <boost/mpi.hpp>
#include <boost/numeric/odeint/util/copy.hpp>
#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include "mpi_state.hpp"

namespace boost {
namespace numeric {
namespace odeint {


/** \brief Copy data from some container on node 0 to the slaves.
 * SourceContainer must support `s::value_type`, `s::const_iterator`, `s.begin()`, `s.end()` and `s.size()`,
 * with Random Access Iterators; i.e. it must be a Random Access Container. */
template< class SourceContainer >
struct copy_impl< SourceContainer, mpi_state< std::vector< typename SourceContainer::value_type > > >
{
    static void copy( const SourceContainer &from, mpi_state< std::vector< typename SourceContainer::value_type > > &to )
    {
        typedef typename SourceContainer::const_iterator it_t;
        typedef typename SourceContainer::value_type T;
        std::vector< std::vector<T> > pieces;
        if(to.world.rank() == 0) {
            // split as evenly as possible; length difference is at most 1.
            const size_t num = static_cast<size_t>(to.world.size());
            const size_t part = from.size() / num;
            const size_t mod = from.size() % num;
            it_t begin = from.begin(), end = begin;
            pieces.resize(num);
            for(size_t i = 0 ; i < num ; i++) {
                const size_t offset = i < mod ? 1 : 0;
                end = begin + (part + offset);
                pieces[i].resize(end - begin);
                std::copy(begin, end, pieces[i].begin());
                begin = end;
            }
        }
        // send to nodes
        boost::mpi::scatter(to.world, pieces, to.data, 0);
    }
};

/** \brief Copy data from an mpi_state to some container on node 0.
 * TargetContainer must support `s::value_type`, `s::iterator`, `s.begin()` and `s.resize(n)`,
 * i.e. it must be a `std::vector`. */
template< class TargetContainer >
struct copy_impl< mpi_state< std::vector< typename TargetContainer::value_type > >, TargetContainer >
{
    static void copy( const mpi_state< std::vector< typename TargetContainer::value_type > > &from , TargetContainer &to )
    {
        typedef typename TargetContainer::value_type T;
        std::vector< std::vector<T> > pieces;
        // send data to root
        boost::mpi::gather(from.world, from.data, pieces, 0);
        if(from.world.rank() == 0) {
            // resize target
            size_t total_size = 0;
            for(size_t i = 0 ; i < pieces.size() ; i++)
                total_size += pieces[i].size();
            to.resize( total_size );
            // copy parts
            typename TargetContainer::iterator out = to.begin();
            for(size_t i = 0 ; i < pieces.size() ; i++)
                out = std::copy(pieces[i].begin(), pieces[i].end(), out);
        }
    }
};



}
}
}


#endif

