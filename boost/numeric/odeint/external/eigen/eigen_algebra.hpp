#ifndef BOOST_NUMERIC_ODEINT_EXTERNAL_EIGEN_EIGEN_THRUST_ALGEBRA_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_EXTERNAL_EIGEN_EIGEN_THRUST_ALGEBRA_HPP_INCLUDED

#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>

// Necessary routines for Eigen matrices to work with vector_space_algebra
//  from odeint
//  (that is, it lets odeint treat the eigen matrices correctly, knowing
//   how to add, multiply, find max element, etc)

namespace Eigen {
template<typename D>
inline const
typename Eigen::CwiseUnaryOp<
          typename Eigen::internal::scalar_add_op<
               typename Eigen::internal::traits<D>::Scalar>,
          const D >
operator+(const typename Eigen::MatrixBase<D> &m,
          const typename Eigen::internal::traits<D>::Scalar &s) {
     return Eigen::CwiseUnaryOp<
          typename Eigen::internal::scalar_add_op<
               typename Eigen::internal::traits<D>::Scalar>,
          const D >(m.derived(),Eigen::internal::scalar_add_op<
                    typename Eigen::internal::traits<D>::Scalar>(s));
}

template<typename D>
inline const
typename Eigen::CwiseUnaryOp<
          typename Eigen::internal::scalar_add_op<
               typename Eigen::internal::traits<D>::Scalar>,
          const D >
operator+(const typename Eigen::internal::traits<D>::Scalar &s,
		const typename Eigen::MatrixBase<D> &m) {
     return Eigen::CwiseUnaryOp<
          typename Eigen::internal::scalar_add_op<
               typename Eigen::internal::traits<D>::Scalar>,
          const D >(m.derived(),Eigen::internal::scalar_add_op<
                    typename Eigen::internal::traits<D>::Scalar>(s));
}

template<typename D1,typename D2>
inline const
typename Eigen::CwiseBinaryOp<
		typename Eigen::internal::scalar_quotient_op<
			typename Eigen::internal::traits<D1>::Scalar>,
				const D1, const D2>
operator/(const Eigen::MatrixBase<D1> &x1, const Eigen::MatrixBase<D2> &x2) {
	return x1.cwiseQuotient(x2);
}

/*
 // doesn't work (too general... doesn't get selected and the abs in
 // Eigen goes into an infinite recursion)
template<typename D>
inline const
typename Eigen::CwiseUnaryOp<
	typename Eigen::internal::scalar_add_op<typename Eigen::internal::traits<D>::Scalar>,
	const D>
abs(const Eigen::MatrixBase<D> &m) {
	return m.cwiseAbs();
}
*/

template<typename B,int S1,int S2,int O, int M1, int M2>
Eigen::Matrix<B,S1,S2,O,M1,M2> abs(const Eigen::Matrix<B,S1,S2,O,M1,M2> &m) {
	return m.cwiseAbs();
}

} // end Eigen namespace

namespace boost { namespace numeric { namespace odeint {

template<typename B,int S1,int S2,int O, int M1, int M2>
struct vector_space_reduce<Eigen::Matrix<B,S1,S2,O,M1,M2> > {
	template<class Value, class Op>
	Value operator()(const Eigen::Matrix<B,S1,S2,O,M1,M2> &x,
					Op op, Value init) {
		for(int i=0;i<x.rows();i++) for(int j=0;j<x.cols();j++)
			init = op(init,x(i,j));
		return init;
	}
};

}}} // end boost::numeric::odeint namespace

#endif
