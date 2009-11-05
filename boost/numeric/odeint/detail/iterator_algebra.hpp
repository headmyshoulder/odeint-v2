/* Boost odeint/detail/accumulators.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner
 
 Some algebraic operations for iterators

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
namespace it_algebra { // iterator algebra

    /* computes x_n += alpha * y_n for all n
       where x and y are given by the iterators x_start, x_end, y_start
       make sure that x and y have the same range
     */
    template< class Iterator1 ,
	      class Iterator2 ,
	      class T >
    void scale_and_add( Iterator1 x_start ,
			Iterator1 x_end ,
			Iterator2 y_start ,
			T alpha )
    {
	while( x_start != x_end )
	    (*x_start++) += alpha * (*y_start++);
    }


    /* computes out_n = x_n - y_n for all n
       where x,y, and out are given as iterators
       make sure that x, y and out have the same range
    */
    template< class InputIterator1 ,
	      class InputIterator2 ,
	      class OutputIterator >
    void substract_and_assign( InputIterator1 x_start ,
			       InputIterator1 x_end ,
			       InputIterator2 y_start ,
			       OutputIterator out_start )
    {
	while( x_start != x_end )
	    ( *out_start++ ) = ( *x_start++ ) - ( *y_start++ );
    }


    /* computes out_n = x_n + alpha * y_n for all n
       where x,y and out are given as iterators.
       make sure, that x,y and out have the same range
    */
    template< class InputIterator1 ,
	      class InputIterator2 ,
	      class OutputIterator ,
	      class T >
    void scale_and_add_and_assign( InputIterator1 x_start ,
				   InputIterator1 x_end ,
				   InputIterator2 y_start ,
				   OutputIterator out_start ,
				   T alpha )
    {
	while( x_start != x_end )
	    ( *out_start++) = ( *x_start++ ) + alpha * ( *y_start++ );
    }




    /* computes out_n = alpha2* ( x_n + y_n + alpha1*z_n )for all n
       where x,y,z and out are given as iterators.
       make sure, that x,y,z and out have the same range
    */
    template< class InputIterator1 ,
	      class InputIterator2 ,
	      class InputIterator3 ,
	      class OutputIterator ,
	      class T >
    void scale_and_add_and_add_and_assign( InputIterator1 x_start ,
					   InputIterator1 x_end ,
					   InputIterator2 y_start ,
					   InputIterator3 z_start ,
					   OutputIterator out_start ,
					   T alpha1 ,
					   T alpha2 )
    {
	while( x_start != x_end )
	    (*out_start++) += alpha2 * ( (*x_start++) + (*y_start++) + alpha1*(*z_start++) );
    }

}
}
}
}
}


#endif //BOOST_NUMERIC_ODEINT_DETAIL_ACCUMULATORS_HPP
