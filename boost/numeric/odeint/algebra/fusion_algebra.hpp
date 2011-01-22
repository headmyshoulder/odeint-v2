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


#define BOOST_FUSION_UNFUSED_MAX_ARITY 10
#define BOOST_FUSION_UNFUSED_TYPE_MAX_ARITY 10
#define BOOST_FUSION_INVOKE_MAX_ARITY 10
#define BOOST_FUSION_INVOKE_PROCEDURE_MAX_ARITY 10
#define BOOST_FUSION_INVOKE_FUNCTION_OBJECT_MAX_ARITY 10
#include <boost/fusion/container.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/view.hpp>
#include <boost/fusion/functional.hpp>

#define BOOST_FUNCTIONAL_FORWARD_ADAPTER_MAX_ARITY 9
#include <boost/functional/forward_adapter.hpp>



namespace boost {
namespace numeric {
namespace odeint {


struct fusion_algebra
{
	struct for_each1_impl
	{
		template< class S1 , class Op >
		void operator()( S1 &s1 , Op op ) const
		{
			boost::fusion::for_each( s1 , op );
		}
		typedef void result_type;
	};


	struct for_each2_impl
	{
		template< class S1 , class S2 , class Op >
		void operator()( S1 &s1 , S2 &s2 , Op op ) const
		{
			typedef boost::fusion::vector< S1& , S2& > Sequences;
			Sequences sequences( s1 , s2 );
			boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
		}
		typedef void result_type;
	};


	struct for_each3_impl
	{
		template< class S1 , class S2 , class S3 , class Op >
		void operator()( S1 &s1 , S2 &s2 , S3 &s3 , Op op ) const
		{
			typedef boost::fusion::vector< S1& , S2& , S3& > Sequences;
			Sequences sequences( s1 , s2 , s3 );
			boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
		}
		typedef void result_type;
	};

	struct for_each4_impl
	{
		template< class S1 , class S2 , class S3 , class S4 , class Op >
		void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , Op op ) const
		{
			typedef boost::fusion::vector< S1& , S2& , S3& , S4& > Sequences;
			Sequences sequences( s1 , s2 , s3 , s4 );
			boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
		}
		typedef void result_type;
	};


	struct for_each5_impl
	{
		template< class S1 , class S2 , class S3 , class S4 , class S5 , class Op >
		void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , Op op ) const
		{
			typedef boost::fusion::vector< S1& , S2& , S3& , S4& , S5& > Sequences;
			Sequences sequences( s1 , s2 , s3 , s4 , s5 );
			boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
		}
		typedef void result_type;
	};


	struct for_each6_impl
	{
		template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 , class Op >
		void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , Op op ) const
		{
			typedef boost::fusion::vector< S1& , S2& , S3& , S4& , S5& , S6& > Sequences;
			Sequences sequences( s1 , s2 , s3 , s4 , s5 , s6 );
			boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
		}
		typedef void result_type;
	};


	struct for_each7_impl
	{
		template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 , class S7 , class Op >
		void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , Op op ) const
		{
			typedef boost::fusion::vector< S1& , S2& , S3& , S4& , S5& , S6& , S7& > Sequences;
			Sequences sequences( s1 , s2 , s3 , s4 , s5 , s6 , s7 );
			boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
		}
		typedef void result_type;
	};


	struct for_each8_impl
	{
		template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 , class S7 , class S8 , class Op >
		void operator()( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , S8 &s8 , Op op ) const
		{
			typedef boost::fusion::vector< S1& , S2& , S3& , S4& , S5& , S6& , S7& , S8& > Sequences;
			Sequences sequences( s1 , s2 , s3 , s4 , s5 , s6 , s7 , s8 );
			boost::fusion::for_each( boost::fusion::zip_view< Sequences >( sequences ) , boost::fusion::make_fused( op ) );
		}
		typedef void result_type;
	};


	struct reduce_impl
	{
		template< class Value , class S , class Reduction >
		Value operator()( S &s , Reduction red , Value init) const
		{
			return boost::fusion::accumulate( s , init , red );
		}

		template< class T > struct result;
		template< class F , class T1 , class T2 , class T3 >
		struct result< F( T1 , T2 , T3 ) >
		{
			/*
			 * A Simple
			 * typedef T3 type;
			 * would result in warnings
			 */
			typedef typename boost::remove_reference< T3 >::type type;
		};
	};



	typedef boost::forward_adapter< for_each1_impl , 2 > for_each1;
	typedef boost::forward_adapter< for_each2_impl , 3 > for_each2;
	typedef boost::forward_adapter< for_each3_impl , 4 > for_each3;
	typedef boost::forward_adapter< for_each4_impl , 5 > for_each4;
	typedef boost::forward_adapter< for_each5_impl , 6 > for_each5;
	typedef boost::forward_adapter< for_each6_impl , 7 > for_each6;
	typedef boost::forward_adapter< for_each7_impl , 8 > for_each7;
	typedef boost::forward_adapter< for_each8_impl , 9 > for_each8;
	typedef boost::forward_adapter< reduce_impl , 3 > reduce;
};



} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_VECTOR_SPACE_ALGEBRA_HPP_INCLUDED
