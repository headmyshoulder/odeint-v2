/*
 * algebra.hpp
 *
 *  Created on: Nov 7, 2010
 *      Author: karsten
 */

#ifndef ALGEBRA_HPP_
#define ALGEBRA_HPP_

#include <boost/mpl/int.hpp>
#include <boost/mpl/at.hpp>

#include "convert_value.hpp"

namespace mpl = boost::mpl;

template< class state_type , class coef_tye , class index >
struct algebra
{
	static void do_step( state_type &x_tmp , const state_type &x , const state_type *k_vector , double dt )
	{}
};

template< class state_type , class coef_type >
struct algebra< state_type , coef_type , mpl::int_< 0 > >
{
	static void do_step( state_type &x_tmp , const state_type &x , const state_type *k_vector , double dt )
	{
		for( size_t i=0 ; i<x.size() ; ++i )
		{
			x_tmp[i] = x[i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 0 > >::type >::get_value() * k_vector[0][i];
		}
	}
};

template< class state_type , class coef_type >
struct algebra< state_type , coef_type , mpl::int_< 1 > >
{
	static void do_step( state_type &x_tmp , const state_type &x , const state_type *k_vector , double dt )
	{
		for( size_t i=0 ; i<x.size() ; ++i )
		{
			x_tmp[i] = x[i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 0 > >::type >::get_value() * k_vector[0][i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 1 > >::type >::get_value() * k_vector[1][i];
		}
	}
};

template< class state_type , class coef_type >
struct algebra< state_type , coef_type , mpl::int_< 2 > >
{
	static void do_step( state_type &x_tmp , const state_type &x , const state_type *k_vector , double dt )
	{
		for( size_t i=0 ; i<x.size() ; ++i )
		{
			x_tmp[i] = x[i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 0 > >::type >::get_value() * k_vector[0][i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 1 > >::type >::get_value() * k_vector[1][i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 2 > >::type >::get_value() * k_vector[2][i];
		}

	}
};

template< class state_type , class coef_type >
struct algebra< state_type , coef_type , mpl::int_< 3 > >
{
	static void do_step( state_type &x_tmp , const state_type &x , const state_type *k_vector , double dt )
	{
		for( size_t i=0 ; i<x.size() ; ++i )
		{
			x_tmp[i] = x[i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 0 > >::type >::get_value() * k_vector[0][i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 1 > >::type >::get_value() * k_vector[1][i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 2 > >::type >::get_value() * k_vector[2][i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 3 > >::type >::get_value() * k_vector[3][i];
		}
	}
};

template< class state_type , class coef_type >
struct algebra< state_type , coef_type , mpl::int_< 4 > >
{
	static void do_step( state_type &x_tmp , const state_type &x , const state_type *k_vector , double dt )
	{
		for( size_t i=0 ; i<x.size() ; ++i )
		{
			x_tmp[i] = x[i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 0 > >::type >::get_value() * k_vector[0][i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 1 > >::type >::get_value() * k_vector[1][i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 2 > >::type >::get_value() * k_vector[2][i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 3 > >::type >::get_value() * k_vector[3][i] +
					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 4 > >::type >::get_value() * k_vector[4][i];
		}
	}
};



#endif /* ALGEBRA_HPP_ */
