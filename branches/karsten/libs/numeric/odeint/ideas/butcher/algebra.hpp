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

#include <boost/numeric/odeint/algebra/standard_algebra.hpp>
#include <boost/numeric/odeint/algebra/standard_operations.hpp>

#include "convert_value.hpp"

namespace mpl = boost::mpl;

template< class state_type , class coef_tye , class index >
struct algebra
{
	static void do_step( state_type &x_tmp , const state_type &x , const state_type *k_vector , const double dt )
	{}
};

template< class state_type , class coef_type >
struct algebra< state_type , coef_type , mpl::int_< 0 > >
{
	typedef typename boost::numeric::odeint::standard_algebra< state_type > std_algebra;
	typedef typename boost::numeric::odeint::standard_operations< double > std_op;

	typedef typename state_type::iterator iterator;
	typedef typename state_type::const_iterator const_iterator;

	static void do_step( state_type &x_tmp , const state_type &x , const state_type *k_vector , const double dt )
	{
		const static double a1 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 0 > >::type >::get_value();

		iterator first1 = x_tmp.begin() , last1 = x_tmp.end();
		const_iterator first2 = x.begin() , first3 = k_vector[0].begin();
		while( first1 != last1 )
			*first1++ = *first2++ + a1 * *first3++ ;


//		std_algebra::for_each3( x_tmp , x ,  k_vector[0] ,
//				std_op::scale_sum2( 1.0 , a1 ) );


//		for( size_t i=0 ; i<x.size() ; ++i )
//		{
//			x_tmp[i] = x[i] + a1 * k_vector[0][i];
//		}
	}
};

template< class state_type , class coef_type >
struct algebra< state_type , coef_type , mpl::int_< 1 > >
{
	typedef typename boost::numeric::odeint::standard_algebra< state_type > std_algebra;
	typedef typename boost::numeric::odeint::standard_operations< double > std_op;

	typedef typename state_type::iterator iterator;
	typedef typename state_type::const_iterator const_iterator;

	static void do_step( state_type &x_tmp , const state_type &x , const state_type *k_vector , const double dt )
	{
		const static double a1 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 0 > >::type >::get_value();
		const static double a2 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 1 > >::type >::get_value();

		iterator first1 = x_tmp.begin() , last1 = x_tmp.end();
		const_iterator first2 = x.begin() , first3 = k_vector[0].begin() , first4 = k_vector[1].begin();
		while( first1 != last1 )
			*first1++ = *first2++ + a1 * *first3++ + a2 * *first4++;

//		for( size_t i=0 ; i<x.size() ; ++i )
//		{
//			x_tmp[i] = x[i] +
//					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 0 > >::type >::get_value() * k_vector[0][i] +
//					dt * convert_value< typename mpl::at< coef_type , mpl::int_< 1 > >::type >::get_value() * k_vector[1][i];
//		}
	}
};

template< class state_type , class coef_type >
struct algebra< state_type , coef_type , mpl::int_< 2 > >
{
	typedef typename boost::numeric::odeint::standard_algebra< state_type > std_algebra;
	typedef typename boost::numeric::odeint::standard_operations< double > std_op;

	typedef typename state_type::iterator iterator;
	typedef typename state_type::const_iterator const_iterator;

	static void do_step( state_type &x_tmp , const state_type &x , const state_type *k_vector , const double dt )
	{
//		const static double a1 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 0 > >::type >::get_value();
//		const static double a2 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 1 > >::type >::get_value();
//		const static double a3 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 2 > >::type >::get_value();
//
//		iterator first1 = x_tmp.begin() , last1 = x_tmp.end();
//		const_iterator first2 = x.begin() , first3 = k_vector[0].begin() , first4 = k_vector[1].begin() , first5 = k_vector[2].begin();
//		while( first1 != last1 )
//			*first1++ = *first2++ + a1 * *first3++ + a2 * *first4++ + a3 * *first5++;

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
	typedef typename boost::numeric::odeint::standard_algebra< state_type > std_algebra;
	typedef typename boost::numeric::odeint::standard_operations< double > std_op;

	typedef typename state_type::iterator iterator;
	typedef typename state_type::const_iterator const_iterator;

	static void do_step( state_type &x_tmp , const state_type &x , const state_type *k_vector , const double dt )
	{
//		const static double a1 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 0 > >::type >::get_value();
//		const static double a2 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 1 > >::type >::get_value();
//		const static double a3 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 2 > >::type >::get_value();
//		const static double a4 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 3 > >::type >::get_value();
//
//		iterator first1 = x_tmp.begin() , last1 = x_tmp.end();
//		const_iterator first2 = x.begin() , first3 = k_vector[0].begin() , first4 = k_vector[1].begin() , first5 = k_vector[2].begin();
//		const_iterator first6 = k_vector[3].begin();
//		while( first1 != last1 )
//			*first1++ = *first2++ + a1 * *first3++ + a2 * *first4++ + a3 * *first5++ + a4 * *first6++;

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
	typedef typename boost::numeric::odeint::standard_algebra< state_type > std_algebra;
	typedef typename boost::numeric::odeint::standard_operations< double > std_op;

	typedef typename state_type::iterator iterator;
	typedef typename state_type::const_iterator const_iterator;

	static void do_step( state_type &x_tmp , const state_type &x , const state_type *k_vector , const double dt )
	{
//		const static double a1 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 0 > >::type >::get_value();
//		const static double a2 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 1 > >::type >::get_value();
//		const static double a3 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 2 > >::type >::get_value();
//		const static double a4 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 3 > >::type >::get_value();
//		const static double a5 = dt * convert_value< typename mpl::at< coef_type , mpl::int_< 4 > >::type >::get_value();
//
//		iterator first1 = x_tmp.begin() , last1 = x_tmp.end();
//		const_iterator first2 = x.begin() , first3 = k_vector[0].begin() , first4 = k_vector[1].begin() , first5 = k_vector[2].begin();
//		const_iterator first6 = k_vector[3].begin() , first7 = k_vector[4].begin();
//		while( first1 != last1 )
//			*first1++ = *first2++ + a1 * *first3++ + a2 * *first4++ + a3 * *first5++ + a4 * *first6++ + a5 * *first7++;

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
