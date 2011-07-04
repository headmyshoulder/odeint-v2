/*
 * convert_values.hpp
 *
 *  Created on: Nov 7, 2010
 *      Author: karsten
 */

#ifndef CONVERT_VALUE_HPP_
#define CONVERT_VALUE_HPP_


#include <boost/ratio.hpp>


template< class T >
struct convert_value
{
	static double get_value( void )
	{
		return double( T::value );
	}
};


template< long N , long D >
struct convert_value< boost::ratio< N , D > >
{
	typedef typename boost::ratio< N , D > number;
	static double get_value( void )
	{
		return double( number::num ) / double( number::den );
	}
};


#endif /* CONVERT_VALUE_HPP_ */
