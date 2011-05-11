/*
 * nr_rk4_lorenz.cpp
 *
 *  Created on: May 11, 2011
 *      Author: mario
 */

#include <boost/array.hpp>

#include "rk_performance_test_case.hpp"

const size_t dim = 3;

typedef boost::array< double , dim > state_type;

struct lorenz
{
    inline void operator()( const state_type &x , state_type &dxdt , const double t ) const
    {
        const double sigma = 10.0;
        const double R = 28.0;
        const double b = 8.0 / 3.0;
        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = x[0]*x[1] - b * x[2];
    }
};

template< class System , typename T , size_t dim >
void rk4_step( const System sys , boost::array< T , dim > &x , const double t , const double dt )
{   // fast rk4 implementation adapted from the book 'Numerical Recipes'
    size_t i;
    const double hh = dt*0.5;
    const double h6 = dt/6.0;
    const double th = t+hh;
    boost::array< T , dim > dydx , dym , dyt , yt;

    sys( x , dydx , t );

    for( i=0 ; i<dim ; i++ )
        yt[i] = x[i] + hh*dydx[i];

    sys( yt , dyt , th );
    for( i=0 ; i<dim ; i++ )
        yt[i] = x[i] + hh*dyt[i];

    sys( yt , dym , th );
    for( i=0 ; i<dim ; i++ ) {
        yt[i] = x[i] + dt*dym[i];
        dym[i] += dyt[i];
    }
    sys( yt , dyt , t+dt );
    for( i=0 ; i<dim ; i++ )
        x[i] += h6*( dydx[i] + dyt[i] + 2.0*dym[i] );
}


class nr_wrapper
{
public:
    void reset_init_cond()
    {
        m_x[0] = 10.0 * rand() / RAND_MAX;
        m_x[1] = 10.0 * rand() / RAND_MAX;
        m_x[2] = 10.0 * rand() / RAND_MAX;
        m_t = 0.0;
    }

    inline void do_step( const double dt )
    {
        rk4_step( lorenz() , m_x , m_t , dt );
    }

    double state( const size_t i ) const
    { return m_x[i]; }

private:
    state_type m_x;
    double m_t;
};



int main()
{
    nr_wrapper stepper;

    run( stepper );
}
