/*
 * range_algebra.hpp
 *
 *  Created on: Feb 12, 2011
 *      Author: karsten
 */

#ifndef RANGE_ALGEBRA_HPP_
#define RANGE_ALGEBRA_HPP_

#include <iostream>
#include <algorithm>

#include <boost/range.hpp>
#include <boost/functional/forward_adapter.hpp>

struct algebra1
{
	struct for_each1_impl
	{
		template< class S1 , class Op >
		void operator()( S1 &s1 , Op op ) const
		{
			std::for_each( boost::begin( s1 ) , boost::end( s1 ) , op	);
			tmp2 = 0.11;
		}
		typedef void result_type;
	};

	struct for_each2_impl
	{
		template< class S1 , class Op >
		void operator()( S1 &s1 , Op op ) const
		{
			std::for_each( boost::begin( s1 ) , boost::end( s1 ) , op	);
			std::cout << tmp2 << std::endl;
		}
		typedef void result_type;
	};


	typedef boost::forward_adapter< for_each1_impl , 2 > for_each1;

	double tmp2;

};


// alternativ

struct algebra2
{
	template< class S1 , class Op >
	void for_each1( S1 &s1 , Op op )
	{
		tmp = 123.0;
		std::for_each( boost::begin( s1 ) , boost::end( s1 ) , op	);
	};

	double tmp;
};



#endif /* RANGE_ALGEBRA_HPP_ */
