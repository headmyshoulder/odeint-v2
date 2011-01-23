/*
 boost header: BOOST_NUMERIC_ODEINT/standard_operations.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_STANDARD_OPERATIONS_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_STANDARD_OPERATIONS_HPP_INCLUDED

#include <algorithm>
#include <cmath>      // for std::max

#include <boost/utility/result_of.hpp>
#include <boost/units/quantity.hpp>

namespace boost {
namespace numeric {
namespace odeint {

/*
 * Conversion of boost::units for use in standard_operations::rel_error
 */
namespace detail
{
		template< class T >
		struct get_value_impl
		{
			static T value( const T &t ) { return t; }
			typedef T result_type;
		};

		template< class Unit , class T >
		struct get_value_impl< boost::units::quantity< Unit , T > >
		{
			static T value( const boost::units::quantity< Unit , T > &t ) { return t.value(); }
			typedef T result_type;
		};

		template< class T >
		typename get_value_impl< T >::result_type get_value( const T &t ) { return get_value_impl< T >::value( t ); }



		template< class T , class V >
		struct set_value_impl
		{
			static void set_value( T &t , const V &v ) { t = v; }
		};

		template< class Unit , class T , class V >
		struct set_value_impl< boost::units::quantity< Unit , T > , V >
		{
			static void set_value( boost::units::quantity< Unit , T > &t , const V &v ) { t = boost::units::quantity< Unit , T >::from_value( v ); }
		};

		template< class T , class V >
		void set_value( T &t , const V &v ) { return set_value_impl< T , V >::set_value( t , v ); }
}



/*
 * Notes:
 *
 * * the results structs are needed in order to work with fusion_algebra
 */
struct standard_operations
{

	template< class Fac1 = double >
	struct scale_sum1
	{
		const Fac1 m_alpha1;

		scale_sum1( const Fac1 &alpha1 ) : m_alpha1( alpha1 ) { }

		template< class T1 , class T2 >
		void operator()( T1 &t1 , const T2 &t2 ) const
		{
			t1 = m_alpha1 * t2;
		}

		typedef void result_type;
	};


	template< class Fac1 = double , class Fac2 = Fac1 >
	struct scale_sum2
	{
		const Fac1 m_alpha1;
		const Fac2 m_alpha2;

		scale_sum2( const Fac1 &alpha1 , const Fac2 &alpha2 ) : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) { }

		template< class T1 , class T2 , class T3 >
		void operator()( T1 &t1 , const T2 &t2 , const T3 &t3) const
		{
			t1 = m_alpha1 * t2 + m_alpha2 * t3;
		}

		typedef void result_type;
	};


	template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 >
	struct scale_sum3
	{
		const Fac1 m_alpha1;
		const Fac2 m_alpha2;
		const Fac3 m_alpha3;

		scale_sum3( const Fac1 &alpha1 , const Fac2 &alpha2 , const Fac3 &alpha3 )
			: m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) { }

		template< class T1 , class T2 , class T3 , class T4 >
		void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 ) const
		{
			t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4;
		}

		typedef void result_type;
	};


	template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 >
	struct scale_sum4
	{
		const Fac1 m_alpha1;
		const Fac2 m_alpha2;
		const Fac3 m_alpha3;
		const Fac4 m_alpha4;

		scale_sum4( const Fac1 &alpha1 , const Fac2 &alpha2 , const Fac3 &alpha3 , const Fac4 &alpha4 )
				: m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) { }

		template< class T1 , class T2 , class T3 , class T4 , class T5 >
		void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5) const
		{
			t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5;
		}

		typedef void result_type;
	};


	template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 >
	struct scale_sum5
	{
		const Fac1 m_alpha1;
		const Fac2 m_alpha2;
		const Fac3 m_alpha3;
		const Fac4 m_alpha4;
		const Fac5 m_alpha5;

		scale_sum5( const Fac1 &alpha1 , const Fac2 &alpha2 , const Fac3 &alpha3 , const Fac4 &alpha4 , const Fac5 &alpha5 )
			: m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) { }

		template< class T1 , class T2 , class T3 , class T4 , class T5 , class T6 >
		void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5 , const T6 &t6) const
		{
			t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5 + m_alpha5 * t6;
		}

		typedef void result_type;
	};


	template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 >
	struct scale_sum6
	{
		const Fac1 m_alpha1;
		const Fac2 m_alpha2;
		const Fac3 m_alpha3;
		const Fac4 m_alpha4;
		const Fac5 m_alpha5;
		const Fac6 m_alpha6;

		scale_sum6( const Fac1 &alpha1 , const Fac2 &alpha2 , const Fac3 &alpha3 , const Fac4 &alpha4 , const Fac5 &alpha5 , const Fac6 &alpha6 )
			: m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ){ }

		template< class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 >
		void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5 , const T6 &t6 ,const T7 &t7) const
		{
			t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5 + m_alpha5 * t6 + m_alpha6 * t7;
		}

		typedef void result_type;
	};


	template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 >
	struct scale_sum7
	{
		const Fac1 m_alpha1;
		const Fac2 m_alpha2;
		const Fac3 m_alpha3;
		const Fac4 m_alpha4;
		const Fac5 m_alpha5;
		const Fac6 m_alpha6;
		const Fac7 m_alpha7;

		scale_sum7( const Fac1 &alpha1 , const Fac2 &alpha2 , const Fac3 &alpha3 , const Fac4 &alpha4 ,
				    const Fac5 &alpha5 , const Fac6 &alpha6 , const Fac7 &alpha7 )
			: m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) { }

		template< class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 >
		void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5 , const T6 &t6 , const T7 &t7 , const T8 &t8 ) const
		{
			t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5 + m_alpha5 * t6 + m_alpha6 * t7 + m_alpha7 * t8;
		}

		typedef void result_type;
	};












	/*
	 * for usage in for_each2
	 *
	 * Works with boost::units by eliminating the unit
	 */
	template< class Fac1 = double >
	struct rel_error
	{
		const Fac1 m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt;

		rel_error( const Fac1 &eps_abs , const Fac1 &eps_rel , const Fac1 &a_x , const Fac1 &a_dxdt )
			: m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) , m_a_x( a_x ) , m_a_dxdt( a_dxdt ) { }


		template< class T1 , class T2 , class T3 >
		void operator()( const T1 &t1 , const T2 &t2 , T3 &t3 ) const
		{
			using std::abs;
			using detail::get_value;
			using detail::set_value;
			set_value( t3 , abs( get_value( t3 ) ) / ( m_eps_abs + m_eps_rel * ( m_a_x * abs( get_value( t1 ) ) + m_a_dxdt * abs( get_value( t2 ) ) ) ) );
		}
	};


	/*
	 * for usage in reduce
	 *
	 * ToDo : check if T1, T2 are units and if so convert them to normal floats
	 */
	template< class Fac1 = double >
	struct maximum
	{
		template< class T1 , class T2 >
		Fac1 operator()( const T1 &t1 , const T2 &t2 ) const
		{
			using std::max;
			using std::abs;
			return max( abs( t1 ) , abs( t2 ) );
		}
	};


};


} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_STANDARD_OPERATIONS_HPP_INCLUDED
