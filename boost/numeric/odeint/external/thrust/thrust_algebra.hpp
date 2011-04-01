/*
 boost header: BOOST_NUMERIC_ODEINT/thrust_algebra.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_THRUST_ALGEBRA_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_THRUST_ALGEBRA_HPP_INCLUDED

#include <thrust/device_vector.h>
#include <thrust/for_each.h>
#include <thrust/iterator/zip_iterator.h>

#include <boost/range.hpp>

namespace boost {
namespace numeric {
namespace odeint {



/*
 * The const versions are needed for boost.range to work, i.e.
 * it allows you to do
 * for_each1( make_pair( vec1.begin() , vec1.begin() + 10 ) , op );
 */

struct thrust_algebra
{
	template< class StateType , class Operation >
	static void for_each1( StateType &s , Operation op )
	{
		thrust::for_each( boost::begin(s) , boost::begin(s) , op );
	}

	template< class StateType1 , class StateType2 , class Operation >
	static void for_each2( StateType1 &s1 , StateType2 &s2 , Operation op )
	{
		thrust::for_each(
				thrust::make_zip_iterator( thrust::make_tuple( boost::begin(s1) ,
                                                               boost::begin(s2) ) ) ,
				thrust::make_zip_iterator( thrust::make_tuple( boost::end(s1) ,
                                                               boost::end(s2) ) ) ,
				op);
	}

	template< class StateType1 , class StateType2 , class StateType3 , class Operation >
	static void for_each3( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , Operation op )
	{
		thrust::for_each(
				thrust::make_zip_iterator( thrust::make_tuple( boost::begin(s1) ,
                                                               boost::begin(s2) ,
                                                               boost::begin(s3) ) ) ,
				thrust::make_zip_iterator( thrust::make_tuple( boost::end(s1) ,
                                                               boost::end(s2) ,
                                                               boost::end(s3) ) ) ,
				op);
	}

	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 ,
                class Operation >
    static void for_each4( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType3 &s4 ,
                                Operation op )
    {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple( boost::begin(s1) ,
                                                               boost::begin(s2) ,
                                                               boost::begin(s3) ,
                                                               boost::begin(s4) ) ) ,
                thrust::make_zip_iterator( thrust::make_tuple( boost::end(s1) ,
                                                               boost::end(s2) ,
                                                               boost::end(s3) ,
                                                               boost::end(s4) ) ) ,
                op);
    }

	template< class StateType1 , class StateType2 , class StateType3 ,
               class StateType4 , class StateType5 ,class Operation >
    static void for_each5( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 ,
                            StateType5 &s5 , Operation op )
    {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple( boost::begin(s1) ,
                                                               boost::begin(s2) ,
                                                               boost::begin(s3) ,
                                                               boost::begin(s4) ,
                                                               boost::begin(s5) ) ) ,
                thrust::make_zip_iterator( thrust::make_tuple( boost::end(s1) ,
                                                               boost::end(s2) ,
                                                               boost::end(s3) ,
                                                               boost::end(s4) ,
                                                               boost::end(s5) ) ) ,
                op);
    }

	template< class StateType1 , class StateType2 , class StateType3 ,
               class StateType4 , class StateType5 , class StateType6 , class Operation >
    static void for_each6( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 ,
                             StateType5 &s5 , StateType6 &s6 , Operation op )
    {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple( boost::begin(s1) ,
                                                               boost::begin(s2) ,
                                                               boost::begin(s3) ,
                                                               boost::begin(s4) ,
                                                               boost::begin(s5) ,
                                                               boost::begin(s6) ) ) ,
                thrust::make_zip_iterator( thrust::make_tuple( boost::end(s1) ,
                                                               boost::end(s2) ,
                                                               boost::end(s3) ,
                                                               boost::end(s4) ,
                                                               boost::end(s5) ,
                                                               boost::end(s6) ) ) ,
                op);
    }

	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 ,
               class StateType5 , class StateType6 , class StateType7 , class Operation >
    static void for_each7( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 ,
                              StateType5 &s5 , StateType6 &s6 , StateType7 &s7 , Operation op )
    {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple( boost::begin(s1) ,
                                                               boost::begin(s2) ,
                                                               boost::begin(s3) ,
                                                               boost::begin(s4) ,
                                                               boost::begin(s5) ,
                                                               boost::begin(s6) ,
                                                               boost::begin(s7) ) ) ,
                thrust::make_zip_iterator( thrust::make_tuple( boost::end(s1) ,
                                                               boost::end(s2) ,
                                                               boost::end(s3) ,
                                                               boost::end(s4) ,
                                                               boost::end(s5) ,
                                                               boost::end(s6) ,
                                                               boost::end(s7) ) ) ,
                op);
    }
};


} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_THRUST_ALGEBRA_HPP_INCLUDED
