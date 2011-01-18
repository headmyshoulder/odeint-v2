/*
 boost header: BOOST_NUMERIC_ODEINT/standard_operations.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_STANDARD_OPERATIONS_UNITS_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_STANDARD_OPERATIONS_UNITS_HPP_INCLUDED

#include <algorithm> // for std::max

namespace boost {
namespace numeric {
namespace odeint {



struct standard_operations_units
{
	template< class Fac1 , class Fac2 >
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
	};

	template< class Fac1 , class Fac2 >
	scale_sum2< Fac1 , Fac2 > make_scale_sum2( const Fac1 &alpha1 , const Fac2 &alpha2 )
	{
		return scale_sum2< Fac1 , Fac2 >( alpha1 , alpha2 );
	}

};


} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_STANDARD_OPERATIONS_UNITS_HPP_INCLUDED
