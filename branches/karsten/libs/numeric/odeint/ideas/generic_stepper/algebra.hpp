/*
 * algebra.hpp
 *
 *  Created on: Nov 13, 2010
 *      Author: mario
 */

#ifndef ALGEBRA_HPP_
#define ALGEBRA_HPP_

#include <vector>

#include <boost/array.hpp>

#include <boost/numeric/odeint/algebra/standard_algebra.hpp>
#include <boost/numeric/odeint/algebra/standard_operations.hpp>


template< typename state_type , size_t n >
struct algebra
{

	typedef typename boost::numeric::odeint::standard_algebra< state_type > std_algebra;
	typedef typename boost::numeric::odeint::standard_operations< double > std_op;


	static void foreach( state_type &x_tmp , const state_type &x , const boost::array< double , n > &a ,
                            const state_type *k_vector , const double dt )
	{}

};

template< typename state_type>
struct algebra< state_type , 1 >
{

	typedef typename boost::numeric::odeint::standard_algebra< state_type > std_algebra;
	typedef typename boost::numeric::odeint::standard_operations< double > std_op;

	static void foreach( state_type &x_tmp , const state_type &x , const boost::array< double , 1 > &a ,
                            const state_type *k_vector , const double dt )
	{
		std_algebra::for_each3( x_tmp , x ,  k_vector[0] , std_op::scale_sum2( 1.0 , a[0]*dt ) );
	}

};


template< typename state_type>
struct algebra< state_type , 2 >
{

	typedef typename boost::numeric::odeint::standard_algebra< state_type > std_algebra;
	typedef typename boost::numeric::odeint::standard_operations< double > std_op;

	static void foreach( state_type &x_tmp , const state_type &x , const boost::array< double , 2 > &a ,
                            const state_type *k_vector , const double dt )
	{
		std_algebra::for_each4( x_tmp , x ,  k_vector[0] , k_vector[1] , std_op::scale_sum3( 1.0 , a[0]*dt , a[1]*dt ) );
	}

};


#endif /* ALGEBRA_HPP_ */
