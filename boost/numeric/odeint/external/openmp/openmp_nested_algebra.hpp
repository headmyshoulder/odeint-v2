/*
 [auto_generated]
 boost/numeric/odeint/external/openmp/openmp_nested_algebra.hpp

 [begin_description]
 Nested parallelized algebra for OpenMP.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_EXTERNAL_OPENMP_OPENMP_NESTED_ALGEBRA_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_EXTERNAL_OPENMP_OPENMP_NESTED_ALGEBRA_HPP_INCLUDED

#include <boost/numeric/odeint/algebra/norm_result_type.hpp>
#include <boost/numeric/odeint/util/n_ary_helper.hpp>

namespace boost {
namespace numeric {
namespace odeint {

/** \brief OpenMP-parallelized algebra, wrapping another, non-parallelized algebra.
 * The NestedState must be a Random Access Container where the elements are sub-states
 * which will be processed in parallel.
 *
 * Requires a NestedState with `s::value_type`, `s[n]` and `s.size()`.
 */
template< class InnerAlgebra >
struct openmp_nested_algebra
{

// FIXME: _Pragma is C++11.
#define BOOST_ODEINT_GEN_BODY(n) \
    const size_t len = s0.size(); \
    _Pragma("omp parallel for schedule(runtime)") \
    for( size_t i = 0 ; i < len ; i++ ) \
        InnerAlgebra::for_each##n( \
            BOOST_PP_ENUM_BINARY_PARAMS(n, s, [i] BOOST_PP_INTERCEPT) , \
            op \
        );
BOOST_ODEINT_GEN_FOR_EACH(BOOST_ODEINT_GEN_BODY)
#undef BOOST_ODEINT_GEN_BODY


    template< class NestedState >
    static typename norm_result_type< typename NestedState::value_type >::type norm_inf( const NestedState &s )
    {
        using std::max;
        using std::abs;
        typedef typename norm_result_type< typename NestedState::value_type >::type result_type;
        result_type init = static_cast< result_type >( 0 );
#       pragma omp parallel for reduction(max: init) schedule(dynamic)
        for( size_t i = 0 ; i < s.size() ; i++ )
            init = max( init , InnerAlgebra::norm_inf( s[i] ) );
        return init;
    }

};


}
}
}

#endif
