/*
 * matrix_vector_adjuster.hpp
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_UTIL_MATRIX_VECTOR_ADJUST_SIZE_HPP_
#define BOOST_NUMERIC_ODEINT_UTIL_MATRIX_VECTOR_ADJUST_SIZE_HPP_


namespace boost {
namespace numeric {
namespace odeint {

struct matrix_vector_adjust_size
{
	template< class Vector , class Matrix >
	static bool adjust_size( const Vector &v , Matrix &m )
	{
		if( ( m.size1() != v.size() ) || ( m.size2() != v.size() ) )
		{
			m.resize( v.size() , v.size() );
			return true;
		}
		return false;
	}
};

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_UTIL_MATRIX_VECTOR_ADJUSTER_HPP_ */
