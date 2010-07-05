/*
 boost header: xyz/vector_space_1d.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef VECTOR_SPACE_1D_HPP_INCLUDED
#define VECTOR_SPACE_1D_HPP_INCLUDED

#include <boost/operators.hpp>

template< class T >
struct vector_space_1d :
    boost::additive1< vector_space_1d< T > ,
    boost::additive2< vector_space_1d< T > , T ,
    boost::multiplicative2< vector_space_1d< T > , T
    > > >
{
	typedef T value_type;

	value_type m_x;

	vector_space_1d( void ) : m_x( 0.0 ) {}

    vector_space_1d& operator+=( const vector_space_1d& p )
	{
    	m_x += p.m_x;
	    return *this;
	}

    vector_space_1d& operator-=( const vector_space_1d& p )
	{
    	m_x -= p.m_x;
	    return *this;
	}

    vector_space_1d& operator+=( const value_type& val )
	{
    	m_x += val;
	    return *this;
	}

    vector_space_1d& operator-=( const value_type& val )
	{
    	m_x -= val;
	    return *this;
	}

    vector_space_1d& operator*=( const value_type &val )
	{
    	m_x *= val;
	    return *this;
	}

    vector_space_1d& operator/=( const value_type &val )
	{
    	m_x /= val;
	    return *this;
	}
};

#endif // VECTOR_SPACE_1D_HPP_INCLUDED
