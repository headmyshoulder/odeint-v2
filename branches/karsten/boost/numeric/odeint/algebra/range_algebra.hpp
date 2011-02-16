/*
 boost header: BOOST_NUMERIC_ODEINT/range_algebra.hpp

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

#include <boost/numeric/odeint/algebra/detail/macros.hpp>
#include <boost/numeric/odeint/algebra/detail/for_each.hpp>
#include <boost/numeric/odeint/algebra/detail/reduce.hpp>

namespace boost {
namespace numeric {
namespace odeint {

struct range_algebra
{
	template< class S1 , class Op >
	static void for_each1( S1 &s1 , Op op )
	{
		detail::for_each1( boost::begin( s1 ) , boost::end( s1 ) ,
				op	);
	}

	template< class S1 , class S2 , class Op >
	static void for_each2( S1 &s1 , S2 &s2 , Op op )
	{
		detail::for_each2( boost::begin( s1 ) , boost::end( s1 ) ,
				boost::begin( s2 ) , op	);
	}

	template< class S1 , class S2 , class S3 , class Op >
	static void for_each3( S1 &s1 , S2 &s2 , S3 &s3 , Op op )
	{
		detail::for_each3( boost::begin( s1 ) , boost::end( s1 ) , boost::begin( s2 ) ,	boost::begin( s3 ) , op	);
	}

	template< class S1 , class S2 , class S3 , class S4 , class Op >
	static void for_each4( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , Op op )
	{
		detail::for_each4( boost::begin( s1 ) , boost::end( s1 ) , boost::begin( s2 ) , boost::begin( s3 ) , boost::begin( s4 ) , op );
	}

	template< class S1 , class S2 , class S3 , class S4 , class S5 , class Op >
	static void for_each5( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , Op op )
	{
		detail::for_each5( boost::begin( s1 ) , boost::end( s1 ) , boost::begin( s2 ) ,	boost::begin( s3 ) , boost::begin( s4 ) , boost::begin( s5 ) , op );
	}

	template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 , class Op >
	static void for_each6( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , Op op )
	{
		detail::for_each6( boost::begin( s1 ) , boost::end( s1 ) , boost::begin( s2 ) , boost::begin( s3 ) , boost::begin( s4 ) , boost::begin( s5 ) , boost::begin( s6 ) ,	op	);
	}

	template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class Op >
	static void for_each7( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , Op op )
	{
		detail::for_each7( boost::begin( s1 ) , boost::end( s1 ) , boost::begin( s2 ) , boost::begin( s3 ) , boost::begin( s4 ) , boost::begin( s5 ) , boost::begin( s6 ) , boost::begin( s7 ) , op	);
	}

	template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class S8 , class Op >
	static void for_each8( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , S8 &s8 , Op op )
	{
		detail::for_each8( boost::begin( s1 ) , boost::end( s1 ) , boost::begin( s2 ) , boost::begin( s3 ) , boost::begin( s4 ) , boost::begin( s5 ) , boost::begin( s6 ) , boost::begin( s7 ) , boost::begin( s8 ) , op	);
	}

	template< class Value , class S , class Red >
	static Value reduce( const S &s , Red red , Value init)
	{
		return detail::reduce( boost::begin( s ) , boost::end( s ) , red , init );
	}


};

} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_STANDARD_ALGEBRA_HPP_INCLUDED
