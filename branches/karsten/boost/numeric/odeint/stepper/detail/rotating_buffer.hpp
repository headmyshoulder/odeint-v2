/*
 * rotating_buffer.hpp
 *
 *  Created on: May 16, 2011
 *      Author: karsten
 */

#ifndef ROTATING_BUFFER_HPP_
#define ROTATING_BUFFER_HPP_

#define BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_ROTATING_BUFFER_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_ROTATING_BUFFER_HPP_

#include <boost/array.hpp>

#include <boost/numeric/odeint/util/construct.hpp>
#include <boost/numeric/odeint/util/destruct.hpp>
#include <boost/numeric/odeint/util/copy.hpp>


namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< class T , size_t N >
class rotating_buffer
{
public:

	typedef T value_type;
	const static size_t dim = N;

	rotating_buffer( void )
	{
		initialize();
	}

	rotating_buffer( const rotating_buffer &rb )
	{
		initialize();
		copy( rb );
	}

	~rotating_buffer( void )
	{
		destroy();
	}

	rotating_buffer& operator=( const rotating_buffer &rb )
	{
		copy( rb );
		return *this;
	}

	size_t size( void ) const
	{
		return dim;
	}

	value_type& operator[]( size_t i )
	{
		return m_data[ get_index( i ) ];
	}

	const value_type& operator[]( size_t i ) const
	{
		return m_data[ get_index( i ) ];
	}

	void rotate( void )
	{
		++m_first;
		if( m_first == dim ) m_first = 0;
	}

protected:

	value_type m_data[N];

private:

	void initialize( void )
	{
		for( size_t i=0 ; i<N ; ++i )
			boost::numeric::odeint::construct( m_data[i] );
		m_first = 0;
	}

	void copy( const rotating_buffer &rb )
	{
		for( size_t i=0 ; i<N ; ++i )
			boost::numeric::odeint::copy( rb[i] , m_data[i] );
	}

	void destroy( void )
	{
		for( size_t i=0 ; i<N ; ++i )
			boost::numeric::odeint::destruct( m_data[i] );
	}

	size_t get_index( size_t i ) const
	{
		return ( ( i + m_first ) % dim );
	}

	size_t m_first;

};


} // detail
} // odeint
} // numeric
} // boost


#endif /* BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_ROTATING_BUFFER_HPP_ */
