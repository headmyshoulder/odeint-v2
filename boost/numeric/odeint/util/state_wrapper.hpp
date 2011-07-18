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
    typedef state_wrapper< V , boost::true_type > state_wrapper_type;
    //typedef typename V::value_type value_type;
    typedef boost::true_type is_resizeable;

    V m_v;

    template< class StateIn >
    bool same_size( const StateIn &x ) const
    {
        //standard implementation relies on boost.range
        return ( boost::size( m_v ) == boost::size( x ) );
    }

    template< class StateIn >
    bool resize( const StateIn &x )
    {
        //standard resizing done like for std::vector
        if( !same_size( x ) )
        {
            m_v.resize( boost::size( x ) );
            return true;
        } else
            return false;
    }
};


template< class V >
struct state_wrapper< V , boost::false_type > // without resizing
{
    typedef state_wrapper< V , boost::false_type > state_wrapper_type;
    //typedef typename V::value_type value_type;
    typedef boost::false_type is_resizeable;

    V m_v;

    //no resize method
};

}
}
}


#endif /* STATE_WRAPPER_HPP_ */
