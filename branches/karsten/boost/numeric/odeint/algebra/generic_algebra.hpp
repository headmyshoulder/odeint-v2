/*
 * generic_algebra.hpp
 *
 *  Created on: May 19, 2011
 *      Author: mario
 */

#ifndef GENERIC_ALGEBRA_HPP_
#define GENERIC_ALGEBRA_HPP_

namespace boost {
namespace numeric {
namespace odeint {

/** TODO: use boost::begin ... **/

struct generic_algebra
{
    template< size_t n , class StateOut , class StateIn , class DerivIn , typename T , class StateIn2 >
    inline static void foreach( StateOut &x_tmp , const StateIn &x ,
            const boost::array< T , n > &a ,
			const DerivIn &dxdt , const StateIn2 k_vector[n] , const T dt )
    {
        for( size_t i=0 ; i<x.size() ; ++i )
        {
            x_tmp[i] = x[i] + a[0]*dt*dxdt[i];
            for( size_t j = 1 ; j<n ; ++j )
                x_tmp[i] += a[j]*dt*k_vector[j-1][i];
        }
    }

    template< size_t n , class StateOut , class StateIn , class DerivIn , typename T >
    inline static void foreach( StateOut &x_tmp ,
                const boost::array< T , n > &a ,
				const DerivIn &dxdt , const StateIn k_vector[n] , const T dt )
    {
        for( size_t i=0 ; i<x.size() ; ++i )
        {
            x_tmp[i] = a[0]*dt*dxdt[i];
            for( size_t j = 1 ; j<n ; ++j )
                x_tmp[i] += a[j]*dt*k_vector[j-1][i];
         }
    }
};

}
}
}

#endif /* GENERIC_ALGEBRA */