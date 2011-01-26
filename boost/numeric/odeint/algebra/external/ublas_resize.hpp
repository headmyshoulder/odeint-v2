/*
 boost header: BOOST_NUMERIC_ODEINT/ublas_resize.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_UBLAS_RESIZE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_UBLAS_RESIZE_HPP_INCLUDED

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <boost/type_traits/integral_constant.hpp> //for true_type and false_type

namespace boost {
namespace numeric {
namespace odeint {

/*
 * specialization for boost::numeric::ublas::vector
 */
template< class T , class A >
struct is_resizeable< boost::numeric::ublas::vector< T , A > >
{
    struct type : public boost::true_type { };
    const static bool value = type::value;
};


/*
 * specialization for boost::numeric::ublas::matrix
 */
template< class T , class L , class A >
struct is_resizeable< boost::numeric::ublas::matrix< T , L , A > >
{
	struct type : public boost::true_type { };
	const static bool value = type::value;
};

template< class T , class L , class A >
struct resize_impl< boost::numeric::ublas::matrix< T , L , A > , boost::numeric::ublas::matrix< T , L , A > >
{
	static void resize( const boost::numeric::ublas::matrix< T , L , A > &x1 , boost::numeric::ublas::matrix< T , L , A > &x2 )
	{
		x2.resize( x1.size1() , x1.size2() );
	}
};

template< class T , class L , class A >
struct same_size_impl< boost::numeric::ublas::matrix< T , L , A > , boost::numeric::ublas::matrix< T , L , A > >
{
	static bool same_size( const boost::numeric::ublas::matrix< T , L , A > &x1 , boost::numeric::ublas::matrix< T , L , A > &x2 )
	{
		return ( ( x1.size1() == x2.size1() ) && ( x1.size2() == x2.size2() ) );
	}
};





} // odeint
} // numeric
} // boost

#endif /* BOOST_NUMERIC_ODEINT_UBLAS_RESIZE_HPP_INCLUDED */
