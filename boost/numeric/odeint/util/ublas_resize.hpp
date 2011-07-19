/*
 [auto_generated]
 boost/numeric/odeint/util/ublas_resize.hpp

 [begin_description]
 Resizing for ublas::vector and ublas::matrix
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_UTIL_UBLAS_RESIZE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_UTIL_UBLAS_RESIZE_HPP_INCLUDED


#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

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


/*
 * specialization for vector-matrix resizing
 */
template< class T , class L , class A >
struct same_size_impl< boost::numeric::ublas::matrix< T , L , A > , boost::numeric::ublas::matrix< T , L , A > >
{
    static bool same_size( const boost::numeric::ublas::matrix< T , L , A > &x1 , const boost::numeric::ublas::matrix< T , L , A > &x2 )
    {
        return ( ( x1.size1() == x2.size1() ) && ( x1.size2() == x2.size2() ) );
    }
};


template< class T_V , class A_V , class T_M , class L_M , class A_M >
struct resize_impl< boost::numeric::ublas::vector< T_V , A_V > , boost::numeric::ublas::matrix< T_M , L_M , A_M > >
{
    static void resize( const boost::numeric::ublas::vector< T_V , A_V > &x1 ,
            boost::numeric::ublas::matrix< T_M , L_M , A_M > &x2 )
    {
        x2.resize( x1.size() , x1.size() );
    }
};

template< class T_V , class A_V , class T_M , class L_M , class A_M >
struct same_size_impl< boost::numeric::ublas::vector< T_V , A_V > , boost::numeric::ublas::matrix< T_M , L_M , A_M > >
{
    static bool same_size( const boost::numeric::ublas::vector< T_V , A_V > &x1 ,
            const boost::numeric::ublas::matrix< T_M , L_M , A_M > &x2 )
    {
        return ( ( x1.size() == x2.size1() ) && ( x1.size() == x2.size2() ) );
    }
};

/*
 * specialization for boost::numeric::ublas::permutation_matrix
 */

template< class T , class A >
struct is_resizeable< boost::numeric::ublas::permutation_matrix< T , A > >
{
    struct type : public boost::true_type { };
    const static bool value = type::value;
};

} // namespace odeint
} // namespace numeric
} // namespace boost




/*
 * preparing ublas::matrix for boost::range, such that ublas::matrix can be used in all steppers with the range algebra
 */

namespace boost
{
template< class T , class L , class A >
struct range_mutable_iterator< boost::numeric::ublas::matrix< T , L , A > >
{
    typedef typename boost::numeric::ublas::matrix< T , L , A >::array_type::iterator type;
};

template< class T , class L , class A >
struct range_const_iterator< boost::numeric::ublas::matrix< T , L , A > >
{
    typedef typename boost::numeric::ublas::matrix< T , L , A >::array_type::const_iterator type;
};

} // namespace boost


namespace boost { namespace numeric { namespace ublas {

template< class T , class L , class A >
inline typename matrix< T , L , A >::array_type::iterator
range_begin( matrix< T , L , A > &x )
{
    return x.data().begin();
}

template< class T , class L , class A >
inline typename matrix< T , L , A >::array_type::const_iterator
range_begin( const matrix< T , L , A > &x )
{
    return x.data().begin();
}

template< class T , class L , class A >
inline typename matrix< T , L , A >::array_type::iterator
range_end( matrix< T , L , A > &x )
{
    return x.data().end();
}

template< class T , class L , class A >
inline typename matrix< T , L , A >::array_type::const_iterator
range_end( const matrix< T , L , A > &x )
{
    return x.data().end();
}

} } } // nampespace boost::numeric::ublas


#endif // BOOST_NUMERIC_ODEINT_UTIL_UBLAS_RESIZE_HPP_INCLUDED
