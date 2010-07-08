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

namespace boost {
namespace numeric {
namespace odeint {


/*
 * have to be changed if thrust device_vector or gsl_vector are used
 */
template< class Time >
struct standard_operations
{
	typedef Time time_type;

	struct increment
	{
		time_type m_dt;

		increment( time_type dt ) : m_dt( dt ) { }

		template< class T1 , class T2 >
		void operator()( T1 &t1 , const T2 &t2 ) const
		{
			t1 += m_dt * t2;
		}
	};
};


} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_STANDARD_OPERATIONS_HPP_INCLUDED
