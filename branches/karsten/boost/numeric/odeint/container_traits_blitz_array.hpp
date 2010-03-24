/*
 boost header: numeric/odeint/blitz_container_traits.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Container traits for blitz::Array.

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_CONTAINER_TRAITS_BLITZ_ARRAY_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_CONTAINER_TRAITS_BLITZ_ARRAY_HPP_INCLUDED

#include <cstdlib>
#include "container_traits.hpp"
#include <blitz/array.h>

#include<iostream>

namespace boost {
namespace numeric {
namespace odeint {

    template< typename T , int n >
    struct container_traits< blitz::Array< T , n > >
    {

	typedef blitz::Array< T , n > container_type;
	typedef T value_type;
	typedef typename container_type::iterator iterator;
	typedef typename container_type::const_iterator const_iterator;



        static void resize( const container_type &x , container_type &dxdt )
        {
            //dxdt.resize( x.shape() );
            dxdt.resizeAndPreserve( x.shape() );
        }
        
        static bool same_size( const container_type &x1 , const container_type &x2 )
        {
            for( int d=0; d<x1.dimensions(); d++ )
                if( x1.extent(d) != x2.extent(d) )
                    return false;
            return true;
        }

	static void adjust_size( const container_type &x1 , container_type &x2 )
        {
	    //if( !same_size( x1 , x2 ) ) resize( x1 , x2 );
            x2.resizeAndPreserve( x1.shape() );
	}

	static iterator begin( container_type &x )
	{
	    return x.begin();
	}

	static const_iterator begin( const container_type &x )
	{
	    return x.begin();
	}

	static iterator end( container_type &x )
	{
	    return x.end();
	}

	static const_iterator end( const container_type &x )
	{
	    return x.end();
	}
    };

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif //BOOST_NUMERIC_ODEINT_CONTAINER_TRAITS_BLITZ_ARRAY_HPP_INCLUDED
