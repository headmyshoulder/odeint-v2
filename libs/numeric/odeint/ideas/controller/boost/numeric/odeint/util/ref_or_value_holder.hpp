/*
 [auto_generated]
 boost/numeric/odeint/util/ref_or_value_holder.hpp

 [begin_description]
 Util class for holding either a value or a reference.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_UTIL_REF_OR_VALUE_HOLDER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_UTIL_REF_OR_VALUE_HOLDER_HPP_INCLUDED


namespace boost {
namespace numeric {
namespace odeint {

/*
 * ToDo : find nice names for the types in ref_or_value_holder
 */

template< class T , bool RefOrVariable > struct ref_or_value_holder ;
// template< class T , bool RefOrVariable > struct cref_or_value_holder ;

template< class T >
struct ref_or_value_holder< T , true >
{
    typedef T value_type;
    typedef T& holder_type;
    typedef T& constructor_type;
    typedef T& reference_type;

    holder_type m_t;
    ref_or_value_holder( constructor_type t ) : m_t( t ) { }
    reference_type get( void ) { return m_t; }
};

template< class T >
struct ref_or_value_holder< T , false >
{
    typedef T value_type;
    typedef T holder_type;
    typedef const T& constructor_type;
    typedef T& reference_type;

    holder_type m_t;
    ref_or_value_holder( constructor_type t ) : m_t( t ) { }
    reference_type get( void ) { return m_t; }
};

//template< class T >
//struct cref_or_value_holder< T , true >
//{
//    const T &m_t;
//    cref_or_value_holder( const T &t ) : m_t( t ) { }
//    const T& get( void ) { return m_t; }
//};
//
//template< class T >
//struct cref_or_value_holder< T , false >
//{
//    const T m_t;
//    cref_or_value_holder( const T &t ) : m_t( t ) { }
//    const T& get( void ) { return m_t; }
//};


} // namespace odeint
} // namespace numeric
} // namespace boost



#endif // BOOST_NUMERIC_ODEINT_UTIL_REF_OR_VALUE_HOLDER_HPP_INCLUDED
