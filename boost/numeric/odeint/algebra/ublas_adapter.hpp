/*
  [auto_generated]
  boost/numeric/odeint/algebra/ublas_adapter.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_ALGEBRA_UBLAS_ADAPTER_HPP_DEFINED
#define BOOST_NUMERIC_ODEINT_ALGEBRA_UBLAS_ADAPTER_HPP_DEFINED

#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

#include <algorithm>

namespace boost {
namespace numeric {
namespace odeint {


template< class T >
struct vector_space_reduce< boost::numeric::ublas::vector< T >  >
{
    template< class Op >
    T operator()( const boost::numeric::ublas::vector< T > &x , Op op , T init ) const
    {
        return std::accumulate( x.begin() , x.end() , init , op );
    }
};


} // namespace odeint
} // namespace numeric
} // namespace boost




namespace boost {
namespace numeric {
namespace ublas {

template<class T>
struct scalar_abs : public scalar_unary_functor< T >
{
    typedef typename scalar_unary_functor< T >::argument_type argument_type;
    typedef typename scalar_unary_functor< T >::result_type result_type;

    static BOOST_UBLAS_INLINE
    result_type apply (argument_type t)
    {
        using std::abs;
        return abs( t );
    }
};


template<class T1, class T2>
struct scalar_max : public scalar_binary_functor<T1, T2>
{
    typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
    typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
    typedef typename scalar_binary_functor<T1, T2>::result_type result_type;

    static BOOST_UBLAS_INLINE
    result_type apply (argument1_type t1, argument2_type t2)
    {
        using std::max;
        return max( t1 , t2 );
    }
};

template<class T1, class T2>
struct scalar_min : public scalar_binary_functor<T1, T2>
{
    typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
    typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
    typedef typename scalar_binary_functor<T1, T2>::result_type result_type;

    static BOOST_UBLAS_INLINE
    result_type apply (argument1_type t1, argument2_type t2)
    {
        using std::min;
        return min( t1 , t2 );
    }
};

template<class E> 
BOOST_UBLAS_INLINE
typename vector_unary_traits<E, scalar_abs<typename E::value_type> >::result_type
abs( const vector_expression<E> &e )
{
    typedef typename vector_unary_traits<E, scalar_abs<typename E::value_type> >::expression_type expression_type;
    return expression_type (e ());
}




template<class E1, class E2>
BOOST_UBLAS_INLINE
typename vector_binary_traits<E1, E2, scalar_min<typename E1::value_type,
                                                 typename E2::value_type> >::result_type
min( const vector_expression<E1> &e1 , const vector_expression<E2> &e2)
{
    typedef typename vector_binary_traits<E1, E2, scalar_min< typename E1::value_type,
                                                              typename E2::value_type> >::expression_type expression_type;
    return expression_type (e1 (), e2 ());
}

template<class E1, class E2>
BOOST_UBLAS_INLINE
typename vector_binary_traits<E1, E2, scalar_max<typename E1::value_type,
                                                 typename E2::value_type> >::result_type
max( const vector_expression<E1> &e1 , const vector_expression<E2> &e2)
{
    typedef typename vector_binary_traits<E1, E2, scalar_max<typename E1::value_type,
                                                             typename E2::value_type> >::expression_type expression_type;
    return expression_type (e1 (), e2 ());
}

template<class E1, class E2>
BOOST_UBLAS_INLINE
typename vector_binary_traits<E1, E2, scalar_divides<typename E1::value_type,
                                                     typename E2::value_type> >::result_type
operator/( const vector_expression<E1> &e1 , const vector_expression<E2> &e2)
{
    typedef typename vector_binary_traits<E1, E2, scalar_divides<typename E1::value_type,
                                                                 typename E2::value_type> >::expression_type expression_type;
    return expression_type (e1 (), e2 ());
}




template<class E1, class T2>
BOOST_UBLAS_INLINE
typename enable_if< is_convertible<T2, typename E1::value_type >,    
                    typename vector_binary_scalar2_traits<E1, const T2, scalar_plus<typename E1::value_type, T2> >::result_type
                    >::type
operator + (const vector_expression<E1> &e1, const T2 &e2) {
        typedef typename vector_binary_scalar2_traits<E1, const T2, scalar_plus<typename E1::value_type, T2> >::expression_type expression_type;
        return expression_type (e1 (), e2);
}




template<class T1, class E2>
BOOST_UBLAS_INLINE
typename enable_if< is_convertible<T1, typename E2::value_type >,    
                    typename vector_binary_scalar1_traits<const T1, E2, scalar_plus<T1, typename E2::value_type> >::result_type
                    >::type
operator + (const T1 &e1,
            const vector_expression<E2> &e2) {
    typedef typename vector_binary_scalar1_traits<const T1, E2, scalar_plus<T1, typename E2::value_type> >::expression_type expression_type;
    return expression_type (e1, e2 ());
}


}
}
}





#endif // BOOST_NUMERIC_ODEINT_ALGEBRA_UBLAS_ADAPTER_HPP_DEFINED
