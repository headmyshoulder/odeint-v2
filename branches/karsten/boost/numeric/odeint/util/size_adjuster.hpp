/*
 boost header: NUMERIC_ODEINT/adjust_size.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_ADJUST_SIZE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_ADJUST_SIZE_HPP_INCLUDED

#include <algorithm>
//#include <iostream>

#include <boost/noncopyable.hpp>
#include <boost/array.hpp>

#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resize.hpp>
//#include <boost/numeric/odeint/util/default_adjust_size.hpp>


namespace boost {
namespace numeric {
namespace odeint {



/*
 * Tags to specify resize behavior of steppers
 */
struct adjust_size_manually_tag {};
struct adjust_size_initially_tag {};
struct adjust_size_always_tag {};





/*
 * Adjust size functionality with policies and resizeability
 */
template< class Container , size_t Dim >
class size_adjuster : boost::noncopyable
{
public:

	typedef Container container_type;
	static const size_t dim = Dim;

	size_adjuster() : m_is_initialized( false ) , m_states()
	{
		std::fill( m_states.begin() , m_states.end() , static_cast< container_type* >( 0 ) );
	}

	template< class State >
	bool adjust_size( const State &x )
	{
		return adjust_size_by_resizeability( x , typename is_resizeable< container_type >::type() );
	}

	template< class State >
	bool adjust_size_by_policy( const State &x , adjust_size_manually_tag )
	{
	    return false;
	}

	template< class State >
	bool adjust_size_by_policy( const State &x , adjust_size_initially_tag )
	{
		if( !m_is_initialized )
		{
			m_is_initialized = true;
			return adjust_size_by_resizeability( x , typename is_resizeable< container_type >::type() );
		}
		else
		{
		    return false;
		}
	}

	template< class State >
	bool adjust_size_by_policy( const State &x , adjust_size_always_tag )
	{
		return adjust_size_by_resizeability( x , typename is_resizeable< container_type >::type() );
	}

	void register_state( size_t idx , container_type &x )
	{
		m_states[idx] = &x;
	}


private:

	template< class State >
	bool adjust_size_by_resizeability( const State &x , boost::true_type )
	{
	    bool resized = false;
		for( size_t i=0 ; i<dim ; ++i )
		{
            resized |= adjust_size_impl( x , *(m_states[i]) );
		}
		return resized;
	}

	template< class State >
	bool adjust_size_by_resizeability( const State &x , boost::false_type )
	{
	    return false;
	}

	/* adjust size implementation - resizes only if sizes aren't equal */
	template< class Container1 , class Container2 >
	static bool adjust_size_impl( const Container1 &x1 , Container2 &x2 )
	{
	    if( !same_size( x1 , x2 ) )
	    {
	        resize( x1 , x2 );
	        return true;
	    }
	    else
	    {
	        return false;
	    }
	}


private :

	bool m_is_initialized;
	boost::array< container_type* , dim > m_states;
};



} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_ADJUST_SIZE_HPP_INCLUDED
