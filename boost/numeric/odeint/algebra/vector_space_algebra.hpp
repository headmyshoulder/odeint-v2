/*
 boost header: BOOST_NUMERIC_ODEINT/vector_space_algebra.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_BOOST_NUMERIC_ODEINT_VECTOR_SPACE_ALGEBRA_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_VECTOR_SPACE_ALGEBRA_HPP_INCLUDED

#include <boost/type_traits/remove_reference.hpp>


namespace boost {
namespace numeric {
namespace odeint {


/*
 * This class template has to be overload in order to call vector_space_algebra::reduce
 */
template< class State > struct vector_space_reduce;

/*
 * Example:
 */
//template< class LorenzState >
//class vector_space_reduce
//{
//	template< class Value , class Op >
//	Value operator()( const LorenzState &s , Op op , Value init ) const
//	{
//		init = op( init , s.x );
//		init = op( init , s.y );
//		init = op( init , s.z );
//		return init;
//	}
//};



struct vector_space_algebra
{
	template< class S1 , class Op >
	static void for_each1( S1 &s1 , Op op )
	{
		// ToDo : build checks, that the +-*/ operators are well defined
		op( s1 );
	}

	template< class S1 , class S2 , class Op >
	static void for_each2( S1 &s1 , S2 &s2 , Op op )
	{
		op( s1 , s2 );
	}

	template< class S1 , class S2 , class S3 , class Op >
	static void for_each3( S1 &s1 , S2 &s2 , S3 &s3 , Op op )
	{
		op( s1 , s2 , s3 );
	}

	template< class S1 , class S2 , class S3 , class S4 , class Op >
	static void for_each4( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , Op op )
	{
		op( s1 , s2 , s3 , s4 );
	}

	template< class S1 , class S2 , class S3 , class S4 , class S5 , class Op >
	static void for_each5( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , Op op )
	{
		op( s1 , s2 , s3 , s4 , s5 );
	}

	template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 , class Op >
	static void for_each6( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , Op op )
	{
		op( s1 , s2 , s3 , s4 , s5 , s6 );
	}

	template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class Op >
	static void for_each7( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , Op op )
	{
		op( s1 , s2 , s3 , s4 , s5 , s6 , s7 );
	}

	template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class S8 , class Op >
	static void for_each8( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , S8 &s8 , Op op )
	{
		op( s1 , s2 , s3 , s4 , s5 , s6 , s7 , s8 );
	}


	template< class Value , class S , class Red >
	static Value reduce( const S &s , Red red , Value init )
	{
		boost::numeric::odeint::vector_space_reduce< S > r;
		return r( s , red , init );
	}
};


} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_VECTOR_SPACE_ALGEBRA_HPP_INCLUDED
