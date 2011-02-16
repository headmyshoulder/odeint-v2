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
#define BOOST_FUNCTIONAL_FORWARD_ADAPTER_MAX_ARITY 9
#include <boost/functional/forward_adapter.hpp>


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
	struct for_each1_impl
	{
		template< class StateType , class Operation >
		void operator()( StateType &s , Operation op ) const
		{
			thrust::for_each( boost::begin(s) , boost::begin(s) , op );
		}
		typedef void result_type;
	};


	struct for_each2_impl
	{

		template< class StateType1 , class StateType2 , class Operation >
		void operator()( StateType1 &s1 , StateType2 &s2 , Operation op ) const
		{
			thrust::for_each(
					thrust::make_zip_iterator( thrust::make_tuple( boost::begin(s1) , boost::begin(s2) ) ) ,
					thrust::make_zip_iterator( thrust::make_tuple( boost::end(s1) , boost::end(s2) ) ) ,
					op);
		}
		typedef void result_type;
	};



	struct for_each3_impl
	{


		template< class StateType1 , class StateType2 , class StateType3 , class Operation >
		void operator()( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , Operation op ) const
		{
			thrust::for_each(
					thrust::make_zip_iterator( thrust::make_tuple( boost::begin(s1) , boost::begin(s2) , boost::begin(s3) ) ) ,
					thrust::make_zip_iterator( thrust::make_tuple( boost::end(s1) , boost::end(s2) , boost::begin(s3) ) ) ,
					op);
		}
		typedef void result_type;
	};

	typedef boost::forward_adapter< for_each1_impl , 2 > for_each1;
	typedef boost::forward_adapter< for_each2_impl , 3 > for_each2;
	typedef boost::forward_adapter< for_each3_impl , 4 > for_each3;
};


} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_THRUST_ALGEBRA_HPP_INCLUDED
