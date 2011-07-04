/*
 * state_wrapper.hpp
 *
 *  Created on: Jul 4, 2011
 *      Author: mario
 */

#ifndef STATE_WRAPPER_HPP_
#define STATE_WRAPPER_HPP_

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/type_traits/integral_constant.hpp>

#include <boost/numeric/odeint/util/is_resizeable.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template< class V , typename resizeable = typename boost::numeric::odeint::is_resizeable< V >::type>
struct state_wrapper;

//two standard implementations, with and without resizing depending on is_resizeable< StateType >

template< class V >
struct state_wrapper< V , boost::true_type > // with resizing
{
    typedef state_wrapper< V , boost::false_type > state_wrapper_type;

    typedef typename V::value_type value_type;

    V m_v;

    state_wrapper() : m_v() { }

    state_wrapper( V v ) : m_v( v )
    { }

    state_wrapper( const state_wrapper_type &x ) : m_v( x.m_v )
    { }

    state_wrapper_type& operator=( state_wrapper_type &x )
    {
        x.m_v = m_v;
        return *this;
    }

    template< class StateIn >
    bool same_size( const StateIn &x ) const
    {
        //standard implementation relies on boost.range
        return ( boost::size( m_v ) == boost::size( x ) );
    }

    template< class StateIn >
    void resize( const StateIn &x )
    {
        //standard resizing done like for std::vector
        m_v.resize( boost::size( x ) );
    }
};


template< class V >
struct state_wrapper< V , boost::false_type > // without resizing
{
    typedef state_wrapper< V , boost::false_type > state_wrapper_type;

    typedef typename V::value_type value_type;

    V m_v;

    state_wrapper() : m_v() { }

    state_wrapper( const V &v ) : m_v( v ) { }

    state_wrapper( const state_wrapper_type &x ) : m_v( x.m_v )
    { }

    //no resize method
};

}
}
}


#endif /* STATE_WRAPPER_HPP_ */
