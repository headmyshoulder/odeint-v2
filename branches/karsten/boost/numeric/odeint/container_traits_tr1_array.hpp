/*
 boost header: numeric/odeint/container_traits_tr1_array.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_CONTAINER_TRAITS_TR1_ARRAY_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_CONTAINER_TRAITS_TR1_ARRAY_HPP_INCLUDED

#include <tr1/array>
#include <boost/numeric/odeint/container_traits.hpp>

namespace boost {
namespace numeric {
namespace odeint {

    // Template Specialization for fixed size array - no resizing can happen
    template< class T , size_t N >
    struct container_traits< std::tr1::array< T , N > >
    {
    public:

	typedef std::tr1::array< T , N > container_type;
	typedef typename container_type::value_type value_type;
	typedef typename container_type::iterator iterator;
	typedef typename container_type::const_iterator const_iterator;


        static void resize( const container_type &x , container_type &dxdt )
        {
            throw; // should never be called
        }

        static const bool same_size(
                const container_type &x1 ,
                const container_type &x2
            )
        {
            return true; // if this was false, the code wouldn't compile
        }

        static void adjust_size(
                const container_type &x1 ,
                container_type &x2
            )
        {
            if( !same_size( x1 , x2 ) ) throw;
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


#endif //BOOST_NUMERIC_ODEINT_CONTAINER_TRAITS_TR1_ARRAY_HPP_INCLUDED
