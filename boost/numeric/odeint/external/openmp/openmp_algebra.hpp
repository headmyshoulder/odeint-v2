/*
 [auto_generated]
 boost/numeric/odeint/external/openmp/openmp_algebra.hpp

 [begin_description]
 Nested parallelized algebra for OpenMP.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_EXTERNAL_OPENMP_OPENMP_ALGEBRA_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_EXTERNAL_OPENMP_OPENMP_ALGEBRA_HPP_INCLUDED

#include <boost/numeric/odeint/algebra/norm_result_type.hpp>
#include <boost/numeric/odeint/util/n_ary_helper.hpp>

namespace boost {
namespace numeric {
namespace odeint {

/** \brief OpenMP-parallelized algebra.
 *
 * Requires `s.size()` and `s[n]`, i.e. a Random Access Container.
 */
struct openmp_algebra
{

// FIXME: _Pragma is C++11.
#define OPENMP_ALGEBRA(n) \
    const size_t len = s0.size(); \
    _Pragma("omp parallel for") \
    for( size_t i = 0 ; i < len ; i++ ) \
        op( BOOST_PP_ENUM_BINARY_PARAMS(n, s, [i] BOOST_PP_INTERCEPT) );
BOOST_ODEINT_GEN_FOR_EACH(OPENMP_ALGEBRA)
#undef OPENMP_ALGEBRA


    template< class S >
    static typename norm_result_type< S >::type norm_inf( const S &s )
    {
        using std::max;
        using std::abs;
        typedef typename norm_result_type< S >::type result_type;
        result_type init = static_cast< result_type >( 0 );
#       pragma omp parallel for reduction(max: init)
        for( size_t i = 0 ; i < s.size() ; ++i )
            init = max( init , abs(s[i]) );
        return init;
    }

};


}
}
}

#endif
