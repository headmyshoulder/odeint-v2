/*
 boost header: NUMERIC_ODEINT/adjust_size.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_ADJUST_SIZE_UNITS_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_ADJUST_SIZE_UNITS_HPP_INCLUDED

#include <algorithm>

#include <boost/noncopyable.hpp>
#include <boost/array.hpp>

#include <boost/numeric/odeint/util/default_adjust_size.hpp>


namespace boost {
namespace numeric {
namespace odeint {








/*
 * Adjust size functionality with policies and resizeability
 */
template< class Deriv , size_t Dim >
class size_adjuster_units : boost::noncopyable
{
public:

	size_adjuster_units() : m_is_initialized( false ) , m_states()
	{
		std::fill( m_states.begin() , m_states.end() , static_cast< Deriv* >( 0 ) );
	}

	template< class State >
	bool adjust_size( const State &x )
	{
		return adjust_size_by_resizeability( x , typename is_resizeable< Deriv >::type() );
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
			return adjust_size_by_resizeability( x , typename is_resizeable< Deriv >::type() );
		}
		else
		{
		    return false;
		}
	}

	template< class State >
	bool adjust_size_by_policy( const State &x , adjust_size_always_tag )
	{
		return adjust_size_by_resizeability( x , typename is_resizeable< Deriv >::type() );
	}

	void register_state( size_t idx , Deriv &x )
	{
		m_states[idx] = &x;
	}


private:

	template< class State >
	bool adjust_size_by_resizeability( const State &x , boost::true_type )
	{
		for( size_t i=0 ; i<Dim ; ++i )
		{
            boost::numeric::odeint::default_adjust_size::adjust_size( x , *(m_states[i]) );
		}
		return ( Dim > 0 );
	}

	template< class State >
	bool adjust_size_by_resizeability( const State &x , boost::false_type )
	{
	    return false;
	}


private :

	bool m_is_initialized;
	boost::array< Deriv* , Dim > m_states;
};


} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_ADJUST_SIZE_UNITS_HPP_INCLUDED
