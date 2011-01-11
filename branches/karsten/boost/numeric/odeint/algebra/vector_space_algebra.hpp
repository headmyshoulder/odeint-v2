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

namespace boost {
namespace numeric {
namespace odeint {


struct vector_space_algebra
{
	template< class StateType1 , class StateType2 , class Operation >
	static void for_each2( StateType1 &s1 , StateType2 &s2 , Operation op )
	{
		// ToDo : build checks, that the +-*/ operators are well defined
		op( s1 , s2 );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class Operation >
	static void for_each3( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , Operation op )
	{
		op( s1 , s2 , s3 );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class Operation >
	static void for_each4( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , Operation op )
	{
		op( s1 , s2 , s3 , s4 );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class StateType5 , class Operation >
	static void for_each5( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , StateType5 &s5 , Operation op )
	{
		op( s1 , s2 , s3 , s4 , s5 );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class StateType5 , class StateType6 , class Operation >
	static void for_each6( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , StateType5 &s5 , StateType6 &s6 , Operation op )
	{
		op( s1 , s2 , s3 , s4 , s5 , s6 );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class StateType5 , class StateType6 ,class StateType7 , class Operation >
	static void for_each7( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , StateType5 &s5 , StateType6 &s6 , StateType7 &s7 , Operation op )
	{
		op( s1 , s2 , s3 , s4 , s5 , s6 , s7 );
	}

	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class StateType5 , class StateType6 ,class StateType7 , class StateType8 , class Operation >
	static void for_each8( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , StateType5 &s5 , StateType6 &s6 , StateType7 &s7 , StateType8 &s8 , Operation op )
	{
		op( s1 , s2 , s3 , s4 , s5 , s6 , s7 , s8 );
	}



	/* ToDo : get ValueType from Container? */
	template< class ValueType , class StateType , class Reduction>
	static ValueType reduce( StateType &s , Reduction red , ValueType init)
	{
		return red( s );
	}
};


} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_VECTOR_SPACE_ALGEBRA_HPP_INCLUDED
