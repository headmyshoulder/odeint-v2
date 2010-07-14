/*
 boost header: xyz/gsl_vector_adaptor.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_GSL_VECTOR_ADAPTOR_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_GSL_VECTOR_ADAPTOR_HPP_INCLUDED

#include <gsl/gsl_vector.h>
#include <boost/range.hpp>

#include <boost/iterator/iterator_facade.hpp>

namespace boost {
namespace numeric {
namespace odeint {


class const_gsl_vector_iterator;

/*
 * defines an iterator for gsl_vector
 */
class gsl_vector_iterator : public boost::iterator_facade< gsl_vector_iterator , double , boost::random_access_traversal_tag >
{
public :

	gsl_vector_iterator( void ): m_p(0) , m_stride( 0 ) { }
	explicit gsl_vector_iterator( gsl_vector &p ) : m_p( p.data ) , m_stride( p.stride ) { }
	gsl_vector_iterator( gsl_vector_iterator &p ) : m_p( p.m_p ) , m_stride( p.m_stride ) { }

private :

	friend class boost::iterator_core_access;
	friend class const_gsl_vector_iterator;

	void increment( void ) { m_p += m_stride; }
	void decrement( void ) { m_p -= m_stride; }
	void advance( ptrdiff_t n ) { m_p += n*m_stride; }
	bool equal( const gsl_vector_iterator &other ) const { return this->m_p == other.m_p; }
	bool equal( const const_gsl_vector_iterator &other ) const;
	double& dereference( void ) const { return *m_p; }

	double *m_p;
	size_t m_stride;
};

/*
 * defines an iterator for gsl_vector
 */
class const_gsl_vector_iterator : public boost::iterator_facade< const_gsl_vector_iterator , double , boost::random_access_traversal_tag >
{
public :

	const_gsl_vector_iterator( void ): m_p(0) , m_stride( 0 ) { }
	explicit const_gsl_vector_iterator( const gsl_vector &p ) : m_p( p.data ) , m_stride( p.stride ) { }
	const_gsl_vector_iterator( const const_gsl_vector_iterator &p ) : m_p( p.m_p ) , m_stride( p.m_stride ) { }
	const_gsl_vector_iterator( const gsl_vector_iterator &p ) : m_p( p.m_p ) , m_stride( p.m_stride ) { }

private :

	friend class boost::iterator_core_access;
	friend class gsl_vector_iterator;

	void increment( void ) { m_p += m_stride; }
	void decrement( void ) { m_p -= m_stride; }
	void advance( ptrdiff_t n ) { m_p += n*m_stride; }
	bool equal( const const_gsl_vector_iterator &other ) const { return this->m_p == other.m_p; }
	bool equal( const gsl_vector_iterator &other ) const { return this->m_p == other.m_p; }
	const double& dereference( void ) const { return *m_p; }

	const double *m_p;
	size_t m_stride;
};

bool gsl_vector_iterator::equal( const const_gsl_vector_iterator &other ) const { return this->m_p == other.m_p; }


} // namespace odeint
} // namespace numeric
} // namespace boost


namespace boost {



/*
 * specialization of range_iterator for gsl_vector_iterator
 */
template<>
struct range_iterator< gsl_vector >
{
	typedef boost::numeric::odeint::gsl_vector_iterator type;
};

/*
 * specialization of range_iterator for const_gsl_vector_iterator
 */
template<>
struct range_iterator< const gsl_vector >
{
	typedef boost::numeric::odeint::const_gsl_vector_iterator type;
};


/*
 * specialization of begin for gsl_vector
 */
template<>
range_iterator< gsl_vector >::type begin( gsl_vector &r )
{
	return boost::numeric::odeint::gsl_vector_iterator( r );
}

/*
 * specialization of begin for const gsl_vector
 */
template<>
range_iterator< const gsl_vector >::type begin( const gsl_vector &r )
{
	return boost::numeric::odeint::const_gsl_vector_iterator( r );
}


/*
 * specialization of end for gsl_vector
 */
template<>
range_iterator< gsl_vector >::type end( gsl_vector &r )
{
	return boost::numeric::odeint::gsl_vector_iterator( r ) + r.size * r.stride ;
}

/*
 * specialization of end for const gsl_vector
 */
template<>
range_iterator< const gsl_vector >::type end( const gsl_vector &r )
{
	return boost::numeric::odeint::const_gsl_vector_iterator( r ) + r.size * r.stride ;
}


} // namespace boost



#endif // BOOST_NUMERIC_ODEINT_GSL_VECTOR_ADAPTOR_HPP_INCLUDED
