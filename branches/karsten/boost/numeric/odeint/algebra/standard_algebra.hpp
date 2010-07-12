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
	static void for_each2( StateType1 &s1 , StateType2 &s2 , Operation op )
	{
		// ToDo : check that number of arguments of the operation is equal 2

		// ToDo : generate macro
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType1 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType2 >::type , container_type >::value ));

		// ToDo : pack into detail namespace
		for_each2(	boost::begin( s1 ) , boost::end( s1 ) ,
					boost::begin( s2 ) , op	);
	}

	// ToDo : pack into namespace detail
	template< class Iterator1 , class Iterator2 , class Operation >
	static void for_each2( Iterator1 first1 , Iterator1 last1 , Iterator2 first2 , Operation op )
	{
		for( ; first1 != last1 ; )
			op( *first1++ , *first2++ );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class Operation >
	static void for_each3( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , Operation op )
	{
		// ToDo : check that number of arguments of the operation is equal 3

		// ToDo : generate macro
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType1 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType2 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType3 >::type , container_type >::value ));

		// ToDo : pack into detail namespace
		for_each3(	boost::begin( s1 ) , boost::end( s1 ) ,
					boost::begin( s2 ) ,
					boost::begin( s3 ) ,
					op	);
	}

	// ToDo : pack into namespace detail
	template< class Iterator1 , class Iterator2 , class Iterator3 , class Operation >
	static void for_each3( Iterator1 first1 , Iterator1 last1 , Iterator2 first2 , Iterator3 first3, Operation op )
	{
		for( ; first1 != last1 ; )
			op( *first1++ , *first2++ , *first3++ );
	}



	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class Operation >
	static void for_each4( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , Operation op )
	{
		// ToDo : check that number of arguments of the operation is equal 3

		// ToDo : generate macro
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType1 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType2 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType3 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType4 >::type , container_type >::value ));

		// ToDo : pack into detail namespace
		for_each4(	boost::begin( s1 ) , boost::end( s1 ) ,
					boost::begin( s2 ) ,
					boost::begin( s3 ) ,
					boost::begin( s4 ) ,
					op	);
	}

	// ToDo : pack into namespace detail
	template< class Iterator1 , class Iterator2 , class Iterator3 , class Iterator4 , class Operation >
	static void for_each4( Iterator1 first1 , Iterator1 last1 , Iterator2 first2 , Iterator3 first3, Iterator4 first4, Operation op )
	{
		for( ; first1 != last1 ; )
			op( *first1++ , *first2++ , *first3++ , *first4++ );
	}



	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class StateType5 , class Operation >
	static void for_each5( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , StateType5 &s5 , Operation op )
	{
		// ToDo : check that number of arguments of the operation is equal 3

		// ToDo : generate macro
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType1 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType2 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType3 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType4 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType5 >::type , container_type >::value ));

		// ToDo : pack into detail namespace
		for_each5(	boost::begin( s1 ) , boost::end( s1 ) ,
					boost::begin( s2 ) ,
					boost::begin( s3 ) ,
					boost::begin( s4 ) ,
					boost::begin( s5 ) ,
					op	);
	}

		// ToDo : pack into namespace detail
	template< class Iterator1 , class Iterator2 , class Iterator3 , class Iterator4 , class Iterator5 , class Operation >
	static void for_each5( Iterator1 first1 , Iterator1 last1 , Iterator2 first2 , Iterator3 first3,
			Iterator4 first4, Iterator5 first5, Operation op )
	{
		for( ; first1 != last1 ; )
			op( *first1++ , *first2++ , *first3++ , *first4++ , *first5++ );
	}



	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class StateType5 , class StateType6 , class Operation >
	static void for_each6( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , StateType5 &s5 , StateType6 &s6 , Operation op )
	{
		// ToDo : check that number of arguments of the operation is equal 6

		// ToDo : generate macro
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType1 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType2 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType3 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType4 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType5 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType6 >::type , container_type >::value ));

		// ToDo : pack into detail namespace
		for_each6(	boost::begin( s1 ) , boost::end( s1 ) ,
					boost::begin( s2 ) ,
					boost::begin( s3 ) ,
					boost::begin( s4 ) ,
					boost::begin( s5 ) ,
					boost::begin( s6 ) ,
					op	);
	}

			// ToDo : pack into namespace detail
	template< class Iterator1 , class Iterator2 , class Iterator3 , class Iterator4 , class Iterator5 , class Iterator6 , class Operation >
	static void for_each6( Iterator1 first1 , Iterator1 last1 , Iterator2 first2 , Iterator3 first3,
				Iterator4 first4, Iterator5 first5, Iterator6 first6 , Operation op )
	{
		for( ; first1 != last1 ; )
			op( *first1++ , *first2++ , *first3++ , *first4++ , *first5++ , *first6++ );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class StateType5 , class StateType6 ,class StateType7 , class Operation >
	static void for_each7( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , StateType5 &s5 , StateType6 &s6 , StateType7 &s7 , Operation op )
	{
		// ToDo : check that number of arguments of the operation is equal 6

		// ToDo : generate macro
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType1 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType2 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType3 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType4 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType5 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType6 >::type , container_type >::value ));
		BOOST_STATIC_ASSERT(( boost::is_same< typename boost::remove_const< StateType7 >::type , container_type >::value ));

			// ToDo : pack into detail namespace
		for_each7(	boost::begin( s1 ) , boost::end( s1 ) ,
					boost::begin( s2 ) ,
					boost::begin( s3 ) ,
					boost::begin( s4 ) ,
					boost::begin( s5 ) ,
					boost::begin( s6 ) ,
					boost::begin( s7 ) ,
					op	);
	}

				// ToDo : pack into namespace detail
	template< class Iterator1 , class Iterator2 , class Iterator3 , class Iterator4 , class Iterator5 , class Iterator6 , class Iterator7 , class Operation >
	static void for_each7( Iterator1 first1 , Iterator1 last1 , Iterator2 first2 , Iterator3 first3,
				Iterator4 first4, Iterator5 first5, Iterator6 first6 , Iterator7 first7 , Operation op )
	{
		for( ; first1 != last1 ; )
			op( *first1++ , *first2++ , *first3++ , *first4++ , *first5++ , *first6++ , *first7++ );
	}

};

} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_STANDARD_ALGEBRA_HPP_INCLUDED
