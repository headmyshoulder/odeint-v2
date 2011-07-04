/*
 * generic_rk_call_algebra.hpp
 *
 *  Created on: May 30, 2011
 *      Author: mario
 */

#ifndef GENERIC_RK_CALL_ALGEBRA_HPP_
#define GENERIC_RK_CALL_ALGEBRA_HPP_

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< size_t StageNumber , class Algebra >
struct generic_rk_call_algebra;

template< class Algebra >
struct generic_rk_call_algebra< 1 , Algebra >
{
    template< class S1 , class S2 , class S3 , class S4 , class Op>
    void operator()( S1 &s1 , S2 &s2 ,  S3 &s3 , S4 s4_array[0] , Op op ) const
    {
        Algebra::for_each3( s1 , s2 , s3 , op );
    }

    template< class S1 , class S2 , class S4 , class Op>
    void operator()( S1 &s1 , S2 &s2 , S4 s4_array[0] , Op op ) const
    {
        Algebra::for_each2( s1 , s2 , op );
    }
};

template< class Algebra >
struct generic_rk_call_algebra< 2 , Algebra >
{
    template< class S1 , class S2 , class S3 , class S4 , class Op>
    void operator()( S1 &s1 , S2 &s2 ,  S3 &s3 , S4 s4_array[1] , Op op ) const
    {
        Algebra::for_each4( s1 , s2 , s3 , s4_array[0] , op );
    }

    template< class S1 , class S2 , class S4 , class Op>
    void operator()( S1 &s1 , S2 &s2 , S4 s4_array[1] , Op op ) const
    {
        Algebra::for_each3( s1 , s2 , s4_array[0] , op );
    }
};


template< class Algebra >
struct generic_rk_call_algebra< 3 , Algebra >
{
    template< class S1 , class S2 , class S3 , class S4 , class Op>
    void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 s4_array[2] , Op op ) const
    {
        Algebra::for_each5( s1 , s2 , s3 , s4_array[0] , s4_array[1] , op );
    }

    template< class S1 , class S2 , class S4 , class Op>
    void operator()( S1 &s1 , S2 &s2 , S4 s4_array[2] , Op op ) const
    {
        Algebra::for_each4( s1 , s2 , s4_array[0] , s4_array[1] , op );
    }
};


template< class Algebra >
struct generic_rk_call_algebra< 4 , Algebra >
{
    template< class S1 , class S2 , class S3 , class S4 , class Op>
    void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 s4_array[3] , Op op ) const
    {
        Algebra::for_each6( s1 , s2 , s3 , s4_array[0] , s4_array[1] , s4_array[2] , op );
    }

    template< class S1 , class S2 , class S4 , class Op>
    void operator()( S1 &s1 , S2 &s2 , S4 s4_array[3] , Op op ) const
    {
        Algebra::for_each5( s1 , s2 , s4_array[0] , s4_array[1] , s4_array[2] , op );
    }
};


template< class Algebra >
struct generic_rk_call_algebra< 5 , Algebra >
{
    template< class S1 , class S2 , class S3 , class S4 , class Op>
    void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 s4_array[4] , Op op ) const
    {
        Algebra::for_each7( s1 , s2 , s3 , s4_array[0] , s4_array[1] , s4_array[2] , s4_array[3] , op );
    }

    template< class S1 , class S2 , class S4 , class Op>
    void operator()( S1 &s1 , S2 &s2 , S4 s4_array[4] , Op op ) const
    {
        Algebra::for_each6( s1 , s2 , s4_array[0] , s4_array[1] , s4_array[2] , s4_array[3] , op );
    }
};

template< class Algebra >
struct generic_rk_call_algebra< 6 , Algebra >
{
    template< class S1 , class S2 , class S3 , class S4 , class Op>
    void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 s4_array[5] , Op op ) const
    {
        Algebra::for_each8( s1 , s2 , s3 , s4_array[0] , s4_array[1] , s4_array[2] , s4_array[3] , s4_array[4] , op );
    }

    template< class S1 , class S2 , class S4 , class Op>
    void operator()( S1 &s1 , S2 &s2 , S4 s4_array[5] , Op op ) const
    {
        Algebra::for_each7( s1 , s2 , s4_array[0] , s4_array[1] , s4_array[2] , s4_array[3] , s4_array[4] , op );
    }
};

}
}
}
}

#endif /* GENERIC_RK_CALL_ALGEBRA_HPP_ */
