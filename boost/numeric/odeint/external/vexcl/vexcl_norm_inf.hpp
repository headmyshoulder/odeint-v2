/*
  [auto_generated]
  boost/numeric/odeint/external/vexcl/vexcl_norm_inf.hpp

  [begin_description]
  vector_space_norm_inf specialization for vexcl
  [end_description]

  Copyright 2009-2013 Karsten Ahnert
  Copyright 2009-2013 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_EXTERNAL_VEXCL_VEXCL_NORM_INF_HPP_DEFINED
#define BOOST_NUMERIC_ODEINT_EXTERNAL_VEXCL_VEXCL_NORM_INF_HPP_DEFINED

#include <map>
#include <algorithm>

#include <vexcl/vector.hpp>
#include <vexcl/multivector.hpp>
#include <vexcl/reductor.hpp>

#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

namespace boost {
namespace numeric {
namespace odeint {

namespace detail {

// norm_inf specializations for vexcl types need an instance of
// vex::Reductor<T> to be available. This function returns reference
// to such an instance:
template <typename T>
const vex::Reductor<T, vex::MAX>&
vexcl_reductor(const std::vector<cl::CommandQueue> &queue)
{
    // We will hold one static reductor per set of queues (or, rather, contexts):
    static std::map< std::vector<cl_context>, vex::Reductor<T, vex::MAX> > cache;

    // Extract OpenCL context handles from command queues:
    std::vector<cl_context> ctx;
    ctx.reserve(queue.size());
    for(auto q = queue.begin(); q != queue.end(); ++q)
        ctx.push_back( vex::qctx(*q)() );

    // See if there is suitable instance of reductor already:
    auto r = cache.find(ctx);

    // If not, create new instance and move it to the cache.
    if (r == cache.end())
        r = cache.insert( std::make_pair(
                    std::move(ctx), vex::Reductor<T, vex::MAX>(queue)
                    ) ).first;

    return r->second;
}

} // namespace detail

// specialization for vexcl vector
template <typename T>
struct vector_space_norm_inf< vex::vector<T> > {
    typedef T result_type;

    T operator()( const vex::vector<T> &x ) const {
        auto max = detail::vexcl_reductor<T>(x.queue_list());

        return max( fabs(x) );
    }
};

// specialization for vexcl multivector
template <typename T, size_t N>
struct vector_space_norm_inf< vex::multivector<T, N> > {
    typedef T result_type;

    T operator()( const vex::multivector<T, N> &x ) const {
        auto max = detail::vexcl_reductor<T>(x.queue_list());

        // Reducing a multivector results in std::array<T, N>:
        auto m = max( fabs(x) );

        // We will need to reduce it even further:
        return *std::max_element(m.begin(), m.end());
    }
};


} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_EXTERNAL_VEXCL_VEXCL_NORM_INF_HPP_DEFINED
