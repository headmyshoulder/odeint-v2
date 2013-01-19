/*
  [auto_generated]
  boost/numeric/odeint/external/eigen/eigen_resize.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2012 Ankur Sinha

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_EXTERNAL_EIGEN_EIGEN_RESIZE_HPP_DEFINED
#define BOOST_NUMERIC_ODEINT_EXTERNAL_EIGEN_EIGEN_RESIZE_HPP_DEFINED


#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resize.hpp>
#include <boost/numeric/odeint/util/same_size.hpp>

#include <Eigen/Dense>

namespace boost {
namespace numeric {
namespace odeint {


template < typename Scalar , int Rows , int Cols , int Options , int MaxRows , int MaxCols >
struct is_resizeable< Eigen::Matrix < Scalar , Rows , Cols , Options , MaxRows , MaxCols > >
{ 
    typedef boost::true_type type;
    const static bool value = type::value;
};


template < typename Scalar , int Rows , int Cols , int Options , int MaxRows , int MaxCols >
struct is_resizeable< Eigen::Array < Scalar , Rows , Cols , Options , MaxRows , MaxCols > >
{ 
    typedef boost::true_type type;
    const static bool value = type::value;
};





template< typename Scalar , int Rows , int Cols , int Options , int MaxRows , int MaxCols >
struct same_size_impl<
    Eigen::Matrix < Scalar , Rows , Cols , Options , MaxRows , MaxCols > ,
    Eigen::Matrix < Scalar , Rows , Cols , Options , MaxRows , MaxCols >  >
{
    static bool same_size( const Eigen::Matrix < Scalar , Rows , Cols , Options , MaxRows , MaxCols > &m1 ,
                           const Eigen::Matrix < Scalar , Rows , Cols , Options , MaxRows , MaxCols > &m2 )
    {
        return ((m1.innerSize () == m2.innerSize ()) && (m1.outerSize() == m2.outerSize()));
    }
};

template< typename Scalar , int Rows , int Cols , int Options , int MaxRows , int MaxCols >
struct same_size_impl<
    Eigen::Array < Scalar ,  Rows ,  Cols ,  Options , MaxRows , MaxCols > ,
    Eigen::Array < Scalar ,  Rows ,  Cols ,  Options , MaxRows , MaxCols > >
{
    static bool same_size( const Eigen::Array < Scalar ,  Rows ,  Cols ,  Options ,  MaxRows ,  MaxCols > &v1  ,
                           const Eigen::Array < Scalar ,  Rows ,  Cols ,  Options ,  MaxRows ,  MaxCols > &v2 )
    {
        return (v1.innerSize () == v2.innerSize ()) && (v1.outerSize() == v2.outerSize());
    }
};




template< typename Scalar , int Rows , int Cols , int Options , int MaxRows , int MaxCols >
struct resize_impl<
    Eigen::Matrix < Scalar ,  Rows ,  Cols ,  Options ,  MaxRows ,  MaxCols > ,
    Eigen::Matrix < Scalar ,  Rows ,  Cols ,  Options ,  MaxRows ,  MaxCols > >
{
    static void resize( Eigen::Matrix < Scalar ,  Rows ,  Cols ,  Options ,  MaxRows ,  MaxCols > &m1 ,
                        Eigen::Matrix < Scalar ,  Rows ,  Cols ,  Options ,  MaxRows ,  MaxCols > &m2 )
    {
        m1.resizeLike(m2);
    }
};

template< typename Scalar , int Rows , int Cols , int Options , int MaxRows , int MaxCols >
struct resize_impl<
    Eigen::Array < Scalar ,  Rows ,  Cols ,  Options ,  MaxRows ,  MaxCols > ,
    Eigen::Array < Scalar ,  Rows ,  Cols ,  Options ,  MaxRows ,  MaxCols > >
{
    static void resize( Eigen::Array < Scalar ,  Rows ,  Cols ,  Options ,  MaxRows ,  MaxCols > &v1  , 
                        Eigen::Array < Scalar ,  Rows ,  Cols ,  Options ,  MaxRows ,  MaxCols > &v2 )
    {
        v1.resizeLike(v2);
    }
};



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_EXTERNAL_EIGEN_EIGEN_RESIZE_HPP_DEFINED
