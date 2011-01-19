/*
 boost header: BOOST_NUMERIC_ODEINT/mkl_operations.hpp

 Algebra for using the Intel Math Kernel Library Blas1 routines in odeint

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_MKL_OPERATIONS_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_MKL_OPERATIONS_HPP_INCLUDED

#include <mkl_blas.h>

/* exemplary example for writing bindings to the Intel MKL library
 * see test/mkl for how to use mkl with odeint
 * this is a quick and dirty implementation showing the general possibility.
 * It works only with containers based on double and sequentiel memory allocation.
 */

namespace boost {
namespace numeric {
namespace odeint {

/* only defined for doubles */
struct mkl_operations
{
	template< class Fac1 , class Fac2 > struct scale_sum2;


	template<>
    struct scale_sum2< double , double >
    {
		typedef double Fac1;
		typedef double Fac2;
		const Fac1 m_alpha1;
        const Fac2 m_alpha2;

        scale_sum2( const Fac1 alpha1 , const Fac2 alpha2 ) : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) { }

        template< class T1 , class T2 , class T3 >
        void operator()( T1 &t1 , const T2 &t2 , const T3 &t3) const
        {
            // we get Containers that have size() and [i]-access
            const int n = t1.size();
            const int one = 1;
            //boost::numeric::odeint::copy( t1 , t3 );
            t1 = t2;
            if ( m_alpha1 != 1.0 )
                dscal( &n , &m_alpha1 , &(t1[0]) , &one );
            daxpy( &n , &m_alpha2 , &(t3[0]) , &one , &(t1[0]) , &one );
        }
    };

};

} // odeint
} // numeric
} // boost

#endif /* BOOST_NUMERIC_ODEINT_MKL_OPERATIONS_HPP_INCLUDED */
