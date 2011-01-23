/*
 boost header: BOOST_NUMERIC_ODEINT/standard_resize.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STANDARD_RESIZE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STANDARD_RESIZE_HPP_INCLUDED

#include <vector>
#include <list>

#include <boost/type_traits/integral_constant.hpp> //for true_type and false_type

namespace boost {
namespace numeric {
namespace odeint {



/*
 * by default any type is not resizable
 */
template< class Container >
struct is_resizeable
{
	struct type : public boost::false_type { };
	const static bool value = type::value;
};

/*
 * specialization for std::vector
 */
template< class V, class A >
struct is_resizeable< std::vector< V , A  > >
{
	struct type : public boost::true_type { };
	const static bool value = type::value;
};


/*
 * specialization for std::list
 */
template< class V , class A >
struct is_resizeable< std::list< V , A > >
{
	struct type : public boost::true_type { };
	const static bool value = type::value;
};



template< class Container >
void construct( Container &x )
{
}

template< class Container >
void destruct( Container &x )
{
}

template< class Container , class Deriv >
void resize( const Container &x , Deriv &dxdt )
{
	dxdt.resize( x.size() );
}

template< class Container , class Deriv >
bool same_size( const Container &x1 , const Deriv &x2 )
{
	return ( x1.size() == x2.size() );
}

template< class Container , class Deriv >
bool adjust_size( const Container &x1 , Deriv &x2 )
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

template< class Container , class Deriv >
void copy( const Container &from , Deriv &to )
{
	to = from;
}



} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_STANDARD_RESIZE_HPP_INCLUDED
