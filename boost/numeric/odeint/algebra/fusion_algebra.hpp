/*
 boost header: BOOST_NUMERIC_ODEINT/vector_space_algebra.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_FUSION_ALGEBRA_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_FUSION_ALGEBRA_HPP_INCLUDED

#include <boost/fusion/container.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/view.hpp>
#include <boost/fusion/functional.hpp>


namespace boost {
namespace numeric {
namespace odeint {


/*
 * TODO :
 * 1. testing, include int unit test
 * 2. change the standard operations, using boost::result_of. for example:
 *
 * struct increment
 * {
 *	template< class T > struct result;
 *
 *	template< class F , class T1 , class T2 >
 *	struct result< F( T1 , T2 ) >
 *	{
 *		 typedef void type;
 *	};
 *
 *	template< class T1 , class T2 >
 *	void operator()( T1 &t1 , T2 &t2 ) const
 *	{
 *		t1 += t2;
 *	}
 * };
 *
 *
 */
struct fusion_algebra
{

	template< class StateType1 , class Operation >
	static void for_each1( StateType1 &s1 , Operation op )
	{
		boost::fusion::for_each( s1 , op );
	}


	template< class StateType1 , class StateType2 , class Operation >
	static void for_each2( StateType1 &s1 , StateType2 &s2 , Operation op )
	{
		typedef boost::fusion::vector< StateType1& , StateType2& > Sequences;
		Sequences sequences( s1 , s2 );
		boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class Operation >
	static void for_each3( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , Operation op )
	{
		typedef boost::fusion::vector< StateType1& , StateType2& , StateType3& > Sequences;
		Sequences sequences( s1 , s2 , s3 );
		boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class Operation >
	static void for_each4( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , Operation op )
	{
		typedef boost::fusion::vector< StateType1& , StateType2& , StateType3& , StateType4& > Sequences;
		Sequences sequences( s1 , s2 , s3 , s4 );
		boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class StateType5 , class Operation >
	static void for_each5( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , StateType5 &s5 , Operation op )
	{
		typedef boost::fusion::vector< StateType1& , StateType2& , StateType3& , StateType4& , StateType5& > Sequences;
		Sequences sequences( s1 , s2 , s3 , s4 , s5 );
		boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class StateType5 , class StateType6 , class Operation >
	static void for_each6( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , StateType5 &s5 , StateType6 &s6 , Operation op )
	{
		typedef boost::fusion::vector< StateType1& , StateType2& , StateType3& , StateType4& , StateType5& , StateType6& > Sequences;
		Sequences sequences( s1 , s2 , s3 , s4 , s5 , s6 );
		boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class StateType5 , class StateType6 , class StateType7 , class Operation >
	static void for_each7( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , StateType5 &s5 , StateType6 &s6 , StateType7 &s7 , Operation op )
	{
		typedef boost::fusion::vector< StateType1& , StateType2& , StateType3& , StateType4& , StateType5& , StateType6& , StateType7& > Sequences;
		Sequences sequences( s1 , s2 , s3 , s4 , s5 , s6 , s7 );
		boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
	}


	template< class StateType1 , class StateType2 , class StateType3 , class StateType4 , class StateType5 , class StateType6 , class StateType7 , class StateType8 , class Operation >
	static void for_each8( StateType1 &s1 , StateType2 &s2 , StateType3 &s3 , StateType4 &s4 , StateType5 &s5 , StateType6 &s6 , StateType7 &s7 , StateType8 &s8 , Operation op )
	{
		typedef boost::fusion::vector< StateType1& , StateType2& , StateType3& , StateType4& , StateType5& , StateType6& , StateType7& , StateType8& > Sequences;
		Sequences sequences( s1 , s2 , s3 , s4 , s5 , s6 , s7 , s8 );
		boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
	}


	template< class ValueType , class StateType , class Reduction >
	static ValueType reduce( StateType &s , Reduction red , ValueType init)
	{
		return boost::fusion::accumulate( s , init , red );
	}
};



} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_VECTOR_SPACE_ALGEBRA_HPP_INCLUDED
