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

#include <boost/range.hpp>

#include <boost/type_traits/integral_constant.hpp> //for true_type and false_type

namespace boost {
namespace numeric {
namespace odeint {


/*
 * Default implementation for constructing a container does nothing
 * gsl_vector must be construct explicitly
 */
template< class Container >
struct construct_impl
{
	static void construct( Container &x )
	{
	}
};


template< class Container >
void construct( Container &x )
{
	construct_impl< Container >::construct( x );
}


/*
 * Default implementation for destruction of a container does nothing
 * gsl_vector must be destroyed explicitly
 */
template< class Container >
struct destruct_impl
{
	static void destruct( Container &x )
	{
	}
};

template< class Container >
void destruct( Container &x )
{
	destruct_impl< Container >::destruct( x );
}



/*
 * Default implementation of the copy operation used the assign operator
 * gsl_vector must copied differently
 */
template< class Container1, class Container2 >
struct copy_impl
{
	static void copy( const Container1 &from , Container2 &to )
	{
		to = from;
	}
};

template< class Container1 , class Container2 >
void copy( const Container1 &from , Container2 &to )
{
	copy_impl< Container1 , Container2 >::copy( from , to );
}







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






/*
 * Default implementation of resize functionality.
 * This struct has to be specialized in order to work with ublas::matrix, etc.
 */
template< class Container1 , class Container2 >
struct resize_impl
{
	static void resize( const Container1 &x1 , Container2 &x2 )
	{
		x2.resize( boost::size( x1 ) );
	}
};

template< class Container1 , class Container2 >
void resize( const Container1 &x1 , Container2 &x2 )
{
	resize_impl< Container1 , Container2 >::resize( x1 , x2 );
}




/*
 * Default implementation of same_size functionality.
 * This struct has to be specialized in order to work with ublas::matrix, etc.
 */
template< class Container1 , class Container2 >
struct same_size_impl
{
	static bool same_size( const Container1 &x1 , const Container2 &x2 )
	{
		return ( boost::size( x1 ) == boost::size( x2 ) );
	}
};

template< class Container1 , class Container2 >
bool same_size( const Container1 &x1 , const Container2 &x2 )
{
	return same_size_impl< Container1 , Container2 >::same_size( x1 , x2 );
}





/*
 * Default implementation of adjust size functionality.
 * This struct can be specialized.
 *
 * Return true or false if the container has been resized.
 */
template< class Container1 , class Container2 >
struct adjust_size_impl
{
	static bool adjust_size( const Container1 &x1 , Container2 &x2 )
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
};


template< class Container1 , class Container2 >
bool adjust_size( const Container1 &x1 , Container2 &x2 )
{
	return adjust_size_impl< Container1 , Container2 >::adjust_size( x1 , x2 );
}




} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_STANDARD_RESIZE_HPP_INCLUDED
