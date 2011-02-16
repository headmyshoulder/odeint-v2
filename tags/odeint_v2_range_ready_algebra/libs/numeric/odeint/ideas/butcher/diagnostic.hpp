/*
 * diagnostic.hpp
 *
 *  Created on: Nov 7, 2010
 *      Author: karsten
 */

#ifndef DIAGNOSTIC_HPP_
#define DIAGNOSTIC_HPP_

#include <iostream>

#include <boost/mpl/for_each.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/size.hpp>

#include "convert_value.hpp"

namespace mpl = boost::mpl;

struct print_value
{
	template< class T >
	void operator()( T )
	{
		std::cout << convert_value< T >::get_value() << " ";
	}
};

struct print_vector
{
	template< class Vec >
	void operator()( Vec )
	{
		mpl::for_each< Vec >( print_value() );
		std::cout << "\n";
	}
};

struct print_tableau_vector
{
	template< class Vec >
	void operator()( Vec )
	{
		typedef typename mpl::at< Vec , mpl::int_< 0 > >::type index;
		typedef typename mpl::at< Vec , mpl::int_< 1 > >::type c;
		typedef typename mpl::at< Vec , mpl::int_< 2 > >::type coef;

		std::cout << "Size : " << mpl::size< Vec >::value << "\t";
		std::cout << convert_value< index >::get_value() << " ";
		std::cout << convert_value< c >::get_value() << " ";
		mpl::for_each< coef >( print_value() );
		std::cout << "\n";
	}
};


#endif /* DIAGNOSTIC_HPP_ */
