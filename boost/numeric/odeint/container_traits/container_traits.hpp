/* Boost odeint/container_traits.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file includes container_traits functionality for containers

 container_traits

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_CONTAINER_TRAITS_HPP
#define BOOST_NUMERIC_ODEINT_CONTAINER_TRAITS_HPP

namespace boost {
namespace numeric {
namespace odeint {

    template< class Container >
    struct container_traits
    {

	typedef Container container_type;

	typedef typename container_type::value_type value_type;
	typedef typename container_type::iterator iterator;
	typedef typename container_type::const_iterator const_iterator;


        static void resize( const container_type &x , container_type &dxdt )
        {
            dxdt.resize( x.size() );
        }
        
        static bool same_size(
                const container_type &x1 ,
                const container_type &x2
            )
        {
            return (x1.size() == x2.size());
        }

	static void adjust_size(
                const container_type &x1 ,
                container_type &x2
            )
        {
	    if( !same_size( x1 , x2 ) ) resize( x1 , x2 );
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


#endif // BOOST_NUMERIC_ODEINT_CONTAINER_TRAITS_HPP
