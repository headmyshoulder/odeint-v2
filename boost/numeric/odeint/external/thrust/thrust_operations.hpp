/*
 boost header: BOOST_NUMERIC_ODEINT/thrust_operations.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_THRUST_OPERATIONS_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_THRUST_OPERATIONS_HPP_INCLUDED

namespace boost {
namespace numeric {
namespace odeint {

#include <thrust/tuple.h>
#include <thrust/iterator/zip_iterator.h>


struct thrust_operations
{
	template< class Fac1 = double , class Fac2 = Fac1 >
	struct scale_sum2
	{
		const Fac1 m_alpha1;
		const Fac2 m_alpha2;

		scale_sum2( const Fac1 alpha1 , const Fac2 alpha2 )
            : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) { }

		template< class Tuple >
		__host__ __device__
		void operator()( Tuple t ) const
		{
			thrust::get<0>(t) = m_alpha1 * thrust::get<1>(t) + m_alpha2 * thrust::get<2>(t);
		}
	};

	template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac1 >
    struct scale_sum3
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;

        scale_sum3( const Fac1 alpha1 , const Fac2 alpha2 , const Fac3 alpha3 )
            : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) { }

        template< class Tuple >
        __host__ __device__
        void operator()( Tuple t ) const
        {
            thrust::get<0>(t) = m_alpha1 * thrust::get<1>(t) +
                                m_alpha2 * thrust::get<2>(t) +
                                m_alpha3 * thrust::get<3>(t);
        }
    };


	template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac1 , class Fac4 = Fac1 >
    struct scale_sum4
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;

        scale_sum4( const Fac1 alpha1 , const Fac2 alpha2 , const Fac3 alpha3 , const Fac4 alpha4 )
            : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ){ }

        template< class Tuple >
        __host__ __device__
        void operator()( Tuple t ) const
        {
            thrust::get<0>(t) = m_alpha1 * thrust::get<1>(t) +
                                m_alpha2 * thrust::get<2>(t) +
                                m_alpha3 * thrust::get<3>(t) +
                                m_alpha4 * thrust::get<4>(t);
        }
    };


    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac1 ,
               class Fac4 = Fac1 , class Fac5 = Fac1>
    struct scale_sum5
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;

        scale_sum5( const Fac1 alpha1 , const Fac2 alpha2 , const Fac3 alpha3 ,
                     const Fac4 alpha4 , const Fac5 alpha5 )
            : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) ,
              m_alpha4( alpha4 ) , m_alpha5( alpha5 ) { }

        template< class Tuple >
        __host__ __device__
        void operator()( Tuple t ) const
        {
            thrust::get<0>(t) = m_alpha1 * thrust::get<1>(t) +
                                m_alpha2 * thrust::get<2>(t) +
                                m_alpha3 * thrust::get<3>(t) +
                                m_alpha4 * thrust::get<4>(t) +
                                m_alpha5 * thrust::get<5>(t);
        }
    };


    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac1 ,
               class Fac4 = Fac1 , class Fac5 = Fac1 , class Fac6 = Fac1>
    struct scale_sum6
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;

        scale_sum6( const Fac1 alpha1 , const Fac2 alpha2 , const Fac3 alpha3 ,
                     const Fac4 alpha4 , const Fac5 alpha5 , const Fac6 alpha6 )
            : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) ,
              m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) { }

        template< class Tuple >
        __host__ __device__
        void operator()( Tuple t ) const
        {
            thrust::get<0>(t) = m_alpha1 * thrust::get<1>(t) +
                                m_alpha2 * thrust::get<2>(t) +
                                m_alpha3 * thrust::get<3>(t) +
                                m_alpha4 * thrust::get<4>(t) +
                                m_alpha5 * thrust::get<5>(t) +
                                m_alpha6 * thrust::get<6>(t);
        }
    };

};

} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_THRUST_OPERATIONS_HPP_INCLUDED
