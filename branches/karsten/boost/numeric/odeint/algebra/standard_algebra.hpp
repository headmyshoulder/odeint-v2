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
#define BOOST_FUNCTIONAL_FORWARD_ADAPTER_MAX_ARITY 9
#include <boost/functional/forward_adapter.hpp>

#include <boost/numeric/odeint/algebra/detail/macros.hpp>
#include <boost/numeric/odeint/algebra/detail/for_each.hpp>
#include <boost/numeric/odeint/algebra/detail/reduce.hpp>

namespace boost {
namespace numeric {
namespace odeint {

struct standard_algebra
{
	struct for_each1_impl
	{
		template< class S1 , class Op >
		void operator()( S1 &s1 , Op op ) const
		{
			detail::for_each1( boost::begin( s1 ) , boost::end( s1 ) ,
					op	);
		}
		typedef void result_type;
	};

	struct for_each2_impl
	{
		template< class S1 , class S2 , class Op >
		void operator()( S1 &s1 , S2 &s2 , Op op ) const
		{
			detail::for_each2( boost::begin( s1 ) , boost::end( s1 ) ,
					boost::begin( s2 ) , op	);
		}
		typedef void result_type;
	};


	struct for_each3_impl
	{

		template< class S1 , class S2 , class S3 , class Op >
		void operator()( S1 &s1 , S2 &s2 , S3 &s3 , Op op ) const
		{
			detail::for_each3( boost::begin( s1 ) , boost::end( s1 ) , boost::begin( s2 ) ,	boost::begin( s3 ) , op	);
		}
		typedef void result_type;
	};





	struct for_each4_impl
	{

		template< class S1 , class S2 , class S3 , class S4 , class Op >
		void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , Op op ) const
		{
			detail::for_each4( boost::begin( s1 ) , boost::end( s1 ) , boost::begin( s2 ) , boost::begin( s3 ) , boost::begin( s4 ) , op );
		}
		typedef void result_type;
	};




	struct for_each5_impl
	{

		template< class S1 , class S2 , class S3 , class S4 , class S5 , class Op >
		void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , Op op ) const
		{
			detail::for_each5( boost::begin( s1 ) , boost::end( s1 ) , boost::begin( s2 ) ,	boost::begin( s3 ) , boost::begin( s4 ) , boost::begin( s5 ) , op );
		}
		typedef void result_type;
	};






	struct for_each6_impl
	{

		template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 , class Op >
		void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , Op op ) const
		{
			detail::for_each6( boost::begin( s1 ) , boost::end( s1 ) , boost::begin( s2 ) , boost::begin( s3 ) , boost::begin( s4 ) , boost::begin( s5 ) , boost::begin( s6 ) ,	op	);
		}
		typedef void result_type;
	};






	struct for_each7_impl
	{

		template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class Op >
		void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , Op op ) const
		{
			detail::for_each7( boost::begin( s1 ) , boost::end( s1 ) , boost::begin( s2 ) , boost::begin( s3 ) , boost::begin( s4 ) , boost::begin( s5 ) , boost::begin( s6 ) , boost::begin( s7 ) , op	);
		}
		typedef void result_type;
	};






	struct for_each8_impl
	{
		template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class S8 , class Op >
		void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , S8 &s8 , Op op ) const
		{
			detail::for_each8( boost::begin( s1 ) , boost::end( s1 ) , boost::begin( s2 ) , boost::begin( s3 ) , boost::begin( s4 ) , boost::begin( s5 ) , boost::begin( s6 ) , boost::begin( s7 ) , boost::begin( s8 ) , op	);
		}
		typedef void result_type;
	};




	struct reduce
	{
		template< class Value , class S , class Red >
		Value operator()( const S &s , Red red , Value init) const
		{
			return detail::reduce( boost::begin( s ) , boost::end( s ) , red , init );
		}

//		template< class T > struct result;
//
//		template< class F , class T1 , class T2 , class T3 >
//		struct result< F( T1 , T2 , T3 ) >
//		{
//			typedef T3 type;
//		};
	};


	typedef boost::forward_adapter< for_each1_impl , 2 > for_each1;
	typedef boost::forward_adapter< for_each2_impl , 3 > for_each2;
	typedef boost::forward_adapter< for_each3_impl , 4 > for_each3;
	typedef boost::forward_adapter< for_each4_impl , 5 > for_each4;
	typedef boost::forward_adapter< for_each5_impl , 6 > for_each5;
	typedef boost::forward_adapter< for_each6_impl , 7 > for_each6;
	typedef boost::forward_adapter< for_each7_impl , 8 > for_each7;
	typedef boost::forward_adapter< for_each8_impl , 9 > for_each8;
//	typedef boost::forward_adapter< reduce_impl , 3 > reduce;

};

} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_STANDARD_ALGEBRA_HPP_INCLUDED
