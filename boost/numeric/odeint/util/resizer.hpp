/*
 * resizer.hpp
 *
 *  Created on: Jul 4, 2011
 *      Author: mario
 */

#ifndef RESIZER_HPP_
#define RESIZER_HPP_

#include <boost/numeric/odeint/util/is_resizeable.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template< class ResizeState , class State >
bool adjust_size_by_resizeability( ResizeState &x , const State &y , boost::true_type )
{
    return x.resize( y );
}

template< class ResizeState , class State >
bool adjust_size_by_resizeability( ResizeState &x , const State &y , boost::false_type )
{
    return false;
}

struct always_resizer
{
    template< class State , class ResizeFunction >
    bool adjust_size( const State &x , ResizeFunction f )
    {
        return f( x );
    }
};


struct initially_resizer
{

    bool m_initialized;

    initially_resizer() : m_initialized( false )
    { }

    template< class State , class ResizeFunction >
    bool adjust_size( const State &x , ResizeFunction f )
    {
        if( !m_initialized )
        {
            m_initialized = true;
            return f( x );
        } else
            return false;
    }
};


struct never_resizer
{
    template< class State , class ResizeFunction >
    bool adjust_size( const State &x , ResizeFunction f )
    {
        return false;
    }
};


}
}
}


#endif /* RESIZER_HPP_ */
