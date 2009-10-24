/* Boost odeint/resizer.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file includes resizer functionality for containers

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_RESIZER_HPP
#define BOOST_NUMERIC_ODEINT_RESIZER_HPP

#include <tr1/array>

namespace boost {
namespace numeric {
namespace odeint {

    template< class ContainerType > 
    class resizer
    {
    public:
        void resize( const ContainerType &x , ContainerType &dxdt ) const
        {
            dxdt.resize( x.size() );
        }
        
        bool same_size( const ContainerType &x1 , ContainerType &x2 ) const
        {
            return (x1.size() == x2.size());
        }
    };

    template< class T , size_t N >
    class resizer< std::tr1::array< T , N > >
    {
    public:
        void resize( const std::tr1::array<T,N> &x ,
		     std::tr1::array<T,N> &dxdt ) const
        {
            throw; // should never be called
        }

        const bool same_size( const std::tr1::array<T,N> &x1 ,
			      std::tr1::array<T,N> &x2 ) const
        {
            return true; // if this was false, the code wouldn't compile
        }
    };


} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_RESIZER_HPP
