/*
 boost header: BOOST_NUMERIC_ODEINT/array_algebra.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_NUMERIC_ODEINT_ARRAY_ALGEBRA_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_ARRAY_ALGEBRA_HPP_INCLUDED


#include <boost/array.hpp>

namespace boost {
namespace numeric {
namespace odeint {

struct array_algebra
{
    template< typename T , size_t dim , class Op >
    static void for_each1( boost::array< T , dim > &s1 , Op op )
    {
        for( size_t i=0 ; i<dim ; ++i )
            op( s1[i] );
    }


    template< typename T , size_t dim , class Op >
    static void for_each2( boost::array< T , dim > &s1 ,
                              const boost::array< T , dim > &s2 , Op op )
    {
        for( size_t i=0 ; i<dim ; ++i )
            op( s1[i] , s2[i] );
    }

    template< typename T , size_t dim , class Op >
    static void for_each3( boost::array< T , dim > &s1 ,
                             const boost::array< T , dim > &s2 ,
                             const boost::array< T , dim > &s3 , Op op )
    {
        for( size_t i=0 ; i<dim ; ++i )
            op( s1[i] , s2[i] , s3[i] );
    }

    template< typename T , size_t dim , class Op >
    static void for_each4( boost::array< T , dim > &s1 ,
                              const boost::array< T , dim > &s2 ,
                              const boost::array< T , dim > &s3 ,
                              const boost::array< T , dim > &s4 , Op op )
    {
        for( size_t i=0 ; i<dim ; ++i )
            op( s1[i] , s2[i] , s3[i] , s4[i] );
    }

    template< typename T , size_t dim , class Op >
    static void for_each5( boost::array< T , dim > &s1 ,
                              const boost::array< T , dim > &s2 ,
                              const boost::array< T , dim > &s3 ,
                              const boost::array< T , dim > &s4 ,
                              const boost::array< T , dim > &s5 , Op op )
    {
        for( size_t i=0 ; i<dim ; ++i )
            op( s1[i] , s2[i] , s3[i] , s4[i] , s5[i] );
    }

    template< typename T , size_t dim , class Op >
    static void for_each6( boost::array< T , dim > &s1 ,
                              const boost::array< T , dim > &s2 ,
                              const boost::array< T , dim > &s3 ,
                              const boost::array< T , dim > &s4 ,
                              const boost::array< T , dim > &s5 ,
                              const boost::array< T , dim > &s6 , Op op )
    {
        for( size_t i=0 ; i<dim ; ++i )
            op( s1[i] , s2[i] , s3[i] , s4[i] , s5[i] , s6[i] );
    }

    template< typename T , size_t dim , class Op >
    static void for_each7( boost::array< T , dim > &s1 ,
                              const boost::array< T , dim > &s2 ,
                              const boost::array< T , dim > &s3 ,
                              const boost::array< T , dim > &s4 ,
                              const boost::array< T , dim > &s5 ,
                              const boost::array< T , dim > &s6 ,
                              const boost::array< T , dim > &s7 , Op op )
    {
        for( size_t i=0 ; i<dim ; ++i )
            op( s1[i] , s2[i] , s3[i] , s4[i] , s5[i] , s6[i] , s7[i] );
    }
};

}
}
}

#endif
