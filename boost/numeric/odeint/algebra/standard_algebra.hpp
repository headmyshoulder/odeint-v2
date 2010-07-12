/*
 boost header: BOOST_NUMERIC_ODEINT/standard_algebra.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_STANDARD_ALGEBRA_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_STANDARD_ALGEBRA_HPP_INCLUDED

#include <boost/range.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template< class Container >
struct standard_algebra
{
	typedef Container container_type;

	template< class StateType1 , class StateType2 , class Operation >
	static void transform2( StateType1 &s1 , StateType2 &s2 , Operation op )
	{
		// ToDo : check that number of arguments of the operation is equal 2

		// ToDo : generate macro
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType1 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType2 >::type , container_type >::value ));

		// ToDo : pack into detail namespace
		transform2(	boost::begin( s1 ) , boost::end( s1 ) ,
					boost::begin( s2 ) , op	);
	}

	// ToDo : pack into namespace detail
	template< class Iterator1 , class Iterator2 , class Operation >
	static void transform2( Iterator1 first1 , Iterator1 last1 , Iterator2 first2 , Operation op )
	{
		for( ; first1 != last1 ; )
			op( *first1++ , *first2++ );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class Operation >
	static void transform3( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , Operation op )
	{
		// ToDo : check that number of arguments of the operation is equal 3

		// ToDo : generate macro
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType1 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType2 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType3 >::type , container_type >::value ));

		// ToDo : pack into detail namespace
		transform2(	boost::begin( s1 ) , boost::end( s1 ) ,
					boost::begin( s2 ) ,
					boost::begin( s3 ) ,
					op	);
	}

	// ToDo : pack into namespace detail
	template< class Iterator1 , class Iterator2 , class Iterator3 , class Operation >
	static void transform3( Iterator1 first1 , Iterator1 last1 , Iterator2 first2 , Iterator3 first3, Operation op )
	{
		for( ; first1 != last1 ; )
			op( *first1++ , *first2++ , *first3++ );
	}

};

} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_STANDARD_ALGEBRA_HPP_INCLUDED
