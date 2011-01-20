/*
 * vector_space_reduce.hpp
 *
 *  Created on: Jan 19, 2011
 *      Author: karsten
 */

#ifndef VECTOR_SPACE_REDUCE_HPP_
#define VECTOR_SPACE_REDUCE_HPP_

namespace boost {
namespace numeric {
namespace odeint {

template< class State > struct vector_space_reduce;

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

} // odeint
} // numeric
} // boost


#endif /* VECTOR_SPACE_REDUCE_HPP_ */
