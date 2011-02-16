/*
 * fusion_algebra.hpp
 *
 *  Created on: Nov 13, 2010
 *      Author: mario
 */

#ifndef FUSION_ALGEBRA_HPP_
#define FUSION_ALGEBRA_HPP_

#include <vector>

#include <boost/array.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>

template< size_t n >
struct fusion_algebra
{

	typedef boost::numeric::odeint::range_algebra std_algebra;
	typedef boost::numeric::odeint::default_operations std_op;

	template< class state_type >
	inline static void foreach( state_type &x_tmp , const state_type &x , const boost::array< double , n > &a ,
                            const state_type *k_vector , const double dt )
	{}

};

template<>
struct fusion_algebra< 1 >
{

	typedef boost::numeric::odeint::range_algebra std_algebra;
	typedef boost::numeric::odeint::default_operations std_op;

	template< class state_type >
	inline static void foreach( state_type &x_tmp , const state_type &x , const boost::array< double , 1 > &a ,
                            const state_type *k_vector , const double dt )
	{
		std_algebra::for_each3()( x_tmp , x ,  k_vector[0] , std_op::scale_sum2<>( 1.0 , a[0]*dt ) );
	}

};


template<>
struct fusion_algebra< 2 >
{

	typedef boost::numeric::odeint::range_algebra std_algebra;
	typedef boost::numeric::odeint::default_operations std_op;

	template< class state_type >
	inline static void foreach( state_type &x_tmp , const state_type &x , const boost::array< double , 2 > &a ,
                            const state_type *k_vector , const double dt )
	{
		std_algebra::for_each4()( x_tmp , x ,  k_vector[0] , k_vector[1] ,
                                std_op::scale_sum3<>( 1.0 , a[0]*dt , a[1]*dt ) );
	}

};


template<>
struct fusion_algebra< 3 >
{

    typedef boost::numeric::odeint::range_algebra std_algebra;
    typedef boost::numeric::odeint::default_operations std_op;

    template< class state_type >
    inline static void foreach( state_type &x_tmp , const state_type &x , const boost::array< double , 3 > &a ,
                            const state_type *k_vector , const double dt )
    {
        std_algebra::for_each5()( x_tmp , x ,  k_vector[0] , k_vector[1] , k_vector[2] ,
                                std_op::scale_sum4<>( 1.0 , a[0]*dt , a[1]*dt , a[2]*dt) );
    }

};



template<>
struct fusion_algebra< 4 >
{

    typedef boost::numeric::odeint::range_algebra std_algebra;
    typedef boost::numeric::odeint::default_operations std_op;

    template< class state_type >
    inline static void foreach( state_type &x_tmp , const state_type &x , const boost::array< double , 4 > &a ,
                            const state_type *k_vector , const double dt )
    {
        std_algebra::for_each6()( x_tmp , x ,  k_vector[0] , k_vector[1] , k_vector[2] , k_vector[3] ,
                                std_op::scale_sum5<>( 1.0 , a[0]*dt , a[1]*dt , a[2]*dt , a[3]*dt ) );
    }

};


#endif /* FUSION_ALGEBRA_HPP_ */
