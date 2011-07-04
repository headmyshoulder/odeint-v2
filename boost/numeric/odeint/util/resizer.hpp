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


struct always_resizer
{

    template< class Stepper , class State >
    bool adjust_size( Stepper& stepper, const State &x )
    {
        return stepper.adjust_size( x );
    }

};


struct initially_resizer
{
    bool m_initialized;

    initially_resizer(): m_initialized( false )
    { }

    template< class Stepper , class State >
    bool adjust_size( Stepper& stepper, const State &x )
    {
        if( !m_initialized )
        {
            m_initialized = true;
            return stepper.adjust_size( x );
        }
        return false;
    }
};


struct never_resizer
{
    template< class Stepper , class State >
    void adjust_size( Stepper& stepper, const State &x )
    { }
};


}
}
}


#endif /* RESIZER_HPP_ */
