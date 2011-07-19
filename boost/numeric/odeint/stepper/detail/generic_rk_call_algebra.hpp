/*
 [auto_generated]
 boost/numeric/odeint/stepper/detail/generic_rk_call_algebra.hpp

 [begin_description]
 Algebra caller for the generic Runge-Kutta methods.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_GENERIC_RK_CALL_ALGEBRA_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_GENERIC_RK_CALL_ALGEBRA_HPP_INCLUDED


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
    void operator()( Algebra &algebra , S1 &s1 , S2 &s2 ,  S3 &s3 , S4 s4_array[0] , Op op ) const
    {
        algebra.for_each3( s1 , s2 , s3 , op );
    }

    template< class S1 , class S2 , class S4 , class Op>
    void operator()( Algebra &algebra , S1 &s1 , S2 &s2 , S4 s4_array[0] , Op op ) const
    {
        algebra.for_each2( s1 , s2 , op );
    }
};

template< class Algebra >
struct generic_rk_call_algebra< 2 , Algebra >
{
    template< class S1 , class S2 , class S3 , class S4 , class Op>
    void operator()( Algebra &algebra , S1 &s1 , S2 &s2 ,  S3 &s3 , S4 s4_array[1] , Op op ) const
    {
        algebra.for_each4( s1 , s2 , s3 , s4_array[0].m_v , op );
    }

    template< class S1 , class S2 , class S4 , class Op>
    void operator()( Algebra &algebra , S1 &s1 , S2 &s2 , S4 s4_array[1] , Op op ) const
    {
        algebra.for_each3( s1 , s2 , s4_array[0].m_v , op );
    }
};


template< class Algebra >
struct generic_rk_call_algebra< 3 , Algebra >
{
    template< class S1 , class S2 , class S3 , class S4 , class Op>
    void operator()( Algebra &algebra , S1 &s1 , S2 &s2 , S3 &s3 , S4 s4_array[2] , Op op ) const
    {
        algebra.for_each5( s1 , s2 , s3 , s4_array[0].m_v , s4_array[1].m_v , op );
    }

    template< class S1 , class S2 , class S4 , class Op>
    void operator()( Algebra &algebra , S1 &s1 , S2 &s2 , S4 s4_array[2] , Op op ) const
    {
        algebra.for_each4( s1 , s2 , s4_array[0].m_v , s4_array[1].m_v , op );
    }
};


template< class Algebra >
struct generic_rk_call_algebra< 4 , Algebra >
{
    template< class S1 , class S2 , class S3 , class S4 , class Op>
    void operator()( Algebra &algebra , S1 &s1 , S2 &s2 , S3 &s3 , S4 s4_array[3] , Op op ) const
    {
        algebra.for_each6( s1 , s2 , s3 , s4_array[0].m_v , s4_array[1].m_v , s4_array[2].m_v , op );
    }

    template< class S1 , class S2 , class S4 , class Op>
    void operator()( Algebra &algebra , S1 &s1 , S2 &s2 , S4 s4_array[3] , Op op ) const
    {
        algebra.for_each5( s1 , s2 , s4_array[0].m_v , s4_array[1].m_v , s4_array[2].m_v , op );
    }
};


template< class Algebra >
struct generic_rk_call_algebra< 5 , Algebra >
{
    template< class S1 , class S2 , class S3 , class S4 , class Op>
    void operator()( Algebra &algebra , S1 &s1 , S2 &s2 , S3 &s3 , S4 s4_array[4] , Op op ) const
    {
        algebra.for_each7( s1 , s2 , s3 , s4_array[0].m_v , s4_array[1].m_v , s4_array[2].m_v , s4_array[3].m_v , op );
    }

    template< class S1 , class S2 , class S4 , class Op>
    void operator()( Algebra &algebra , S1 &s1 , S2 &s2 , S4 s4_array[4] , Op op ) const
    {
        algebra.for_each6( s1 , s2 , s4_array[0].m_v , s4_array[1].m_v , s4_array[2].m_v , s4_array[3].m_v , op );
    }
};

template< class Algebra >
struct generic_rk_call_algebra< 6 , Algebra >
{
    template< class S1 , class S2 , class S3 , class S4 , class Op>
    void operator()( Algebra &algebra , S1 &s1 , S2 &s2 , S3 &s3 , S4 s4_array[5] , Op op ) const
    {
        algebra.for_each8( s1 , s2 , s3 , s4_array[0].m_v , s4_array[1].m_v , s4_array[2].m_v , s4_array[3].m_v , s4_array[4].m_v , op );
    }

    template< class S1 , class S2 , class S4 , class Op>
    void operator()( Algebra &algebra , S1 &s1 , S2 &s2 , S4 s4_array[5] , Op op ) const
    {
        algebra.for_each7( s1 , s2 , s4_array[0].m_v , s4_array[1].m_v , s4_array[2].m_v , s4_array[3].m_v , s4_array[4].m_v , op );
    }
};

}
}
}
}

#endif // BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_GENERIC_RK_CALL_ALGEBRA_HPP_INCLUDED
