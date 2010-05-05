/*
 boost header: numeric/odeint/ublas_matrix_container_traits.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_UBLAS_MATRIX_CONTAINER_TRAITS_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_UBLAS_MATRIX_CONTAINER_TRAITS_HPP_INCLUDED

#include "container_traits.hpp"
#include <boost/numeric/ublas/matrix.hpp>

namespace boost {
namespace numeric {
namespace odeint {


    template<class T, class L, class A>
    struct container_traits< boost::numeric::ublas::matrix< T , L , A > >
    {

	typedef boost::numeric::ublas::matrix< T , L , A > container_type;
	typedef typename container_type::value_type value_type;
	typedef typename container_type::array_type::iterator iterator;
	typedef typename container_type::array_type::const_iterator const_iterator;



        static void resize( const container_type &x , container_type &dxdt )
        {
            dxdt.resize( x.size1() , x.size2() );
        }
        
        static bool same_size(
                const container_type &x1 ,
                const container_type &x2
            )
        {
            return ( ( x1.size1() == x2.size1() ) && ( x1.size2() == x2.size2() ) );
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
            return x.data().begin();
        }

        static const_iterator begin( const container_type &x )
        {
            return x.data().begin();
        }

        static iterator end( container_type &x )
        {
            return x.data().end();
        }
        
        static const_iterator end( const container_type &x )
        {
            return x.data().end();
        }
    };

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif //BOOST_NUMERIC_ODEINT_UBLAS_MATRIX_CONTAINER_TRAITS_HPP_INCLUDED
