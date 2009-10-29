/* Boost odeint/detail/accumulators.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner
 
 Some accumulators for odeint

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_DETAIL_ACCUMULATORS_HPP
#define BOOST_NUMERIC_ODEINT_DETAIL_ACCUMULATORS_HPP



namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

    template< class InputIterator ,
	      class OutputIterator ,
	      class T >
    void multiply_and_add( InputIterator first1 ,
			   InputIterator last1 ,
			   OutputIterator first2 ,
			   T dt )
    {
	while( first1 != last1 )
	    (*first1++) += dt * (*first2++);
    }

    template< class InputIterator1 ,
	      class InputIterator2 ,
	      class OutputIterator >
    void substract_and_assign( InputIterator1 first1 ,
			       InputIterator1 last1 ,
			       InputIterator2 first2 ,
			       OutputIterator first3 )
    {
	while( first1 != last1 )
	    ( *first3++ ) = ( *first1++ ) - ( *first2++ );
    }

}
}
}
}


#endif //BOOST_NUMERIC_ODEINT_DETAIL_ACCUMULATORS_HPP
