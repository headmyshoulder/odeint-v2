/*
 boost header: NUMERIC_ODEINT_STEPPER/util/ublas_wrapper.hpp

 Copyright 2011 Karsten Ahnert
 Copyright 2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_UBLAS_WRAPPER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_UBLAS_WRAPPER_HPP_INCLUDED

#include <boost/type_traits/integral_constant.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>

namespace boost {
namespace numeric {
namespace odeint {

/*
 * specialization for boost::numeric::ublas::vector
 */
template< class T , class A >
struct is_resizeable< boost::numeric::ublas::vector< T , A > >
{
    typedef boost::true_type type;
    const static bool value = type::value;
};


/*
 * specialization for boost::numeric::ublas::matrix
 */
template< class T , class L , class A >
struct is_resizeable< boost::numeric::ublas::matrix< T , L , A > >
{
    typedef boost::true_type type;
	const static bool value = type::value;
};


/*
 * specialization for boost::numeric::ublas::permutation_matrix
 */
template< class T , class A >
struct is_resizeable< boost::numeric::ublas::permutation_matrix< T , A > >
{
    typedef boost::true_type type;
    const static bool value = type::value;
};



/* specialization for matrizes because we need to provide matrix-vector resizing */
template< class T , class L , class A >
struct state_wrapper< boost::numeric::ublas::matrix< T , L , A > , boost::true_type > // with resizing
{
    typedef boost::numeric::ublas::matrix< T , L , A > state_type;
    typedef state_wrapper< state_type > state_wrapper_type;
    //typedef typename V::value_type value_type;
    typedef boost::true_type is_resizeable;

    state_type m_v;

    template< class T_V , class A_V >
    bool same_size( const boost::numeric::ublas::vector< T_V , A_V > &x )
    {
        return ( ( x.size() == m_v.size1() ) && ( x.size() == m_v.size2() ) );
    }

    template< class T_V , class A_V >
    bool resize( const boost::numeric::ublas::vector< T_V , A_V > &x )
    {
        if( !same_size( x ) )
        {
            m_v.resize( boost::size( x ) , boost::size( x ) );
            return true;
        } else
            return false;
    }
};


/* specialization for permutation matrizes because we need to provide matrix-vector resizing */
template< class T , class A >
struct state_wrapper< boost::numeric::ublas::permutation_matrix< T , A > , boost::true_type > // with resizing
{
    typedef boost::numeric::ublas::permutation_matrix< T , A > state_type;
    typedef state_wrapper< state_type > state_wrapper_type;
    //typedef typename V::value_type value_type;
    typedef boost::true_type is_resizeable;

    state_type m_v;

    state_wrapper() : m_v( 1 )
    { }

    template< class T_V , class A_V >
    bool same_size( const boost::numeric::ublas::vector< T_V , A_V > &x )
    {
        return ( x.size() == m_v.size() );
    }

    template< class T_V , class A_V >
    bool resize( const boost::numeric::ublas::vector< T_V , A_V > &x )
    {
        //standard resizing done like for std::vector
        if( !same_size( x ) )
        {
            m_v.resize( boost::size( x ) );
            return true;
        } else
            return false;
    }
};

}
}
}

#endif
