/*
 boost header: numeric/odeint/gram_schmitt.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_GRAM_SCHMITT_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_GRAM_SCHMITT_HPP_INCLUDED

#include <iterator>
#include <algorithm>
#include <numeric>

namespace boost {
namespace numeric {
namespace odeint {

    template< class Iterator , class T >
    void normalize( Iterator first , Iterator last , T norm )
    {
	while( first != last ) *first++ /= norm;
    }

    template< class Iterator , class T >
    void substract_vector( Iterator first1 , Iterator last1 ,
			   Iterator first2 , T val )
    {
	while( first1 != last1 )
	    *first1++ -= val * *first2++;
    }

    template< class StateType , class LyapType >
    void gram_schmitt( StateType &x , LyapType &lyap ,
		       size_t n , size_t num_of_lyap )
    {
	if( !num_of_lyap ) return;
	if( ptrdiff_t( ( num_of_lyap + 1 ) * n ) != std::distance( x.begin() , x.end() ) )
	    throw std::domain_error( "renormalization() : size of state does not match the number of lyapunov exponents." );

	typedef typename StateType::value_type value_type;
	typedef typename StateType::iterator iterator;

	value_type norm[num_of_lyap] , tmp[num_of_lyap];
	iterator first = x.begin() + n;
	iterator beg1 = first , end1 = first + n ;
	
	std::fill( norm , norm+num_of_lyap , 0.0 );

	// normalize first vector
	norm[0] = sqrt( std::inner_product( beg1 , end1 , beg1 , 0.0 ) );
	normalize( beg1 , end1 , norm[0] );

	beg1 += n ;
	end1 += n;

	for( size_t j=1 ; j<num_of_lyap ; ++j , beg1+=n , end1+=n )
	{
	    for( size_t k=0 ; k<j ; ++k )
		tmp[k] = std::inner_product( beg1 , end1 , first + k*n , 0.0 );

	    for( size_t k=0 ; k<j ; ++k )
		substract_vector( beg1 , end1 , first + k*n , tmp[k] );
 
	    // nromalize j-th vector
	    norm[j] = sqrt( std::inner_product( beg1 , end1 , beg1 , 0.0 ) );
	    normalize( beg1 , end1 , norm[j] );
	}

	for( size_t j=0 ; j<num_of_lyap ; j++ )
	    lyap[j] += log( norm[j] );
    }


} // namespace odeint
} // namespace numeric
} // namespace boost

#endif //BOOST_NUMERIC_ODEINT_GRAM_SCHMITT_HPP_INCLUDED
