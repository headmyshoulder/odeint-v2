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


    // computes y += alpha * x1
    template <
	class InOutIterator ,
	class InputIterator ,
	class T
	>
    void increment(
	InOutIterator first1 ,
	InOutIterator last1 ,
	InputIterator first2 ,
	T alpha
	)
    {
	while( first1 != last1 )
	    (*first1++) += alpha * (*first2++);
    }


    // computes y = x1 - x2
    template <
	class OutputIterator ,
	class InputIterator1 ,
	class InputIterator2
	>
    void assign_diff(
	OutputIterator first1 ,
	OutputIterator last1 ,
	InputIterator1 first2 ,
	InputIterator2 first3 )
    {
	while( first1 != last1 )
	    (*first1++) = (*first2++) - (*first3++);
    }

    // computes y = x1 + alpha * x2
    template <
	class OutputIterator ,
	class InputIterator1 ,
	class InputIterator2 ,
	class T
	>
    void assign_sum(
	OutputIterator first1 ,
	OutputIterator last1 ,
	InputIterator1 first2 ,
	InputIterator2 first3 ,
	T alpha )
    {
	while( first1 != last1 )
	    (*first1++) = (*first2++) + alpha * (*first3++);
    }


    // computes y = alpha1 * ( x1 + x2 + alpha2*x3 )
    template <
	class OutputIterator ,
	class InputIterator1 ,
	class InputIterator2 ,
	class InputIterator3 ,
	class T
	>
    void increment_sum_sum(
	OutputIterator first1 ,
	OutputIterator last1 ,
	InputIterator1 first2 ,
	InputIterator2 first3 ,
	InputIterator3 first4 ,
	T alpha1 ,
	T alpha2
	)
    {
	while( first1 != last1 )
	    (*first1++) += alpha1 *
		( (*first2++) + (*first3++) + alpha2*(*first4++) );
    }


    // computes y = x1 + alpha * x2 ; x2 += x3
    template<
	class OutputIterator ,
	class InputIterator1 ,
	class InOutIterator ,
	class InputIterator2 ,
	class T
	>
    void assign_sum_increment(
	OutputIterator first1 ,
	OutputIterator last1 ,
	InputIterator1 first2 ,
	InOutIterator first3 ,
	InputIterator2 first4 ,
	T alpha
	)
    {
	while( first1 != last1 )
	{
	    (*first1++) = (*first2++) + alpha * (*first3);
	    (*first3++) += (*first4++);
	}
    }


    




}
}
}
}
}


#endif //BOOST_NUMERIC_ODEINT_DETAIL_ACCUMULATORS_HPP
