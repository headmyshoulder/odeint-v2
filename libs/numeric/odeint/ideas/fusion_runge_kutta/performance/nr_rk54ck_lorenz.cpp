/*
 * nr_rk54ck_lorenz.cpp
 *
 *  Created on: May 12, 2011
 *      Author: mario
 */

#include <boost/array.hpp>

#include "rk_performance_test_case.hpp"

#include "lorenz.hpp"

const size_t dim = 3;

typedef boost::array< double , dim > state_type;


template< class System , typename T , size_t dim >
void rk54ck_step( System sys , boost::array< T , dim > &x , const double t , const double dt , boost::array< T , dim > &xerr )
{   // fast rk54ck implementation adapted from the book 'Numerical Recipes'
    size_t i;
    static const double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875 ,
            b21 = 0.2 , b31 = 3.0 / 40.0 , b32 = 9.0 / 40.0 , b41 = 0.3 , b42 =
                    -0.9 , b43 = 1.2 , b51 = -11.0 / 54.0 , b52 = 2.5 , b53 =
                    -70.0 / 27.0 , b54 = 35.0 / 27.0 , b61 = 1631.0 / 55296.0 ,
            b62 = 175.0 / 512.0 , b63 = 575.0 / 13824.0 , b64 = 44275.0
                    / 110592.0 , b65 = 253.0 / 4096.0 , c1 = 37.0 / 378.0 , c3 =
                    250.0 / 621.0 , c4 = 125.0 / 594.0 , c6 = 512.0 / 1771.0 ,
            dc1 = c1 - 2825.0 / 27648.0 , dc3 = c3 - 18575.0 / 48384.0 , dc4 =
                    c4 - 13525.0 / 55296.0 , dc5 = -277.00 / 14336.0 , dc6 = c6
                    - 0.25;
    const size_t n = dim;
    boost::array< T , dim > dydx , ak2 , ak3 , ak4 , ak5 , ak6 , ytemp;

    sys( x , dydx , t );
    for (i=0;i<n;i++)
        ytemp[i] = x[i] + b21 * dt * dydx[i];

    sys( ytemp , ak2 , t+a2*dt );
    for (i=0;i<n;i++)
        ytemp[i] = x[i] + dt*(b31*dydx[i]+b32*ak2[i]);

    sys( ytemp , ak3 , t+a3*dt );
    for (i=0;i<n;i++)
        ytemp[i] = x[i] + dt*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);

    sys( ytemp , ak4 , t+a4*dt );
    for (i=0;i<n;i++)
        ytemp[i] = x[i] + dt*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);

    sys( ytemp, ak5 , t+a5*dt );
    for (i=0;i<n;i++)
        ytemp[i] = x[i] + dt*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);

    sys( ytemp , ak6 , t+a6*dt );
    for (i=0;i<n;i++)
        x[i] += dt*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
    for (i=0;i<n;i++)
        xerr[i] = dt*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
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
        rk54ck_step( lorenz() , m_x , m_t , dt , m_x_err );
    }

    double state( const size_t i ) const
    { return m_x[i]; }

private:
    state_type m_x;
    state_type m_x_err;
    double m_t;
};



int main()
{
    nr_wrapper stepper;

    run( stepper );
}
