/*
 * lorenz.hpp
 *
 *  Created on: May 12, 2011
 *      Author: mario
 */

#ifndef LORENZ_HPP_
#define LORENZ_HPP_

struct lorenz
{
    template< class state_type >
    void inline operator()( const state_type &x , state_type &dxdt , const double t ) const
    {
        const double sigma = 10.0;
        const double R = 28.0;
        const double b = 8.0 / 3.0;
        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = x[0]*x[1] - b * x[2];
    }
};


template< class state_type >
inline void lorenz_func( const state_type &x , state_type &dxdt , const double t )
{
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

#include <gsl/gsl_matrix.h>

int lorenz_gsl( const double t , const double y[] , double f[] , void *params)
{
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;

    f[0] = sigma * ( y[1] - y[0] );
    f[1] = R * y[0] - y[1] - y[0] * y[2];
    f[2] = y[0]*y[1] - b * y[2];
    return GSL_SUCCESS;
}


#endif /* LORENZ_HPP_ */
