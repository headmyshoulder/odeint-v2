/*
 * gsl_rk4_lorenz.cpp
 *
 *  Created on: May 11, 2011
 *      Author: mario
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "rk_performance_test_case.hpp"

const size_t dim = 3;

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

class gsl_wrapper
{
public:

    gsl_wrapper()
    {
        m_s = gsl_odeiv_step_alloc( gsl_odeiv_step_rk4 , dim);
        m_sys.function = lorenz_gsl;
        m_sys.jacobian = 0;
        m_sys.dimension = dim;
        m_sys.params = 0;
    }

    void reset_init_cond()
    {
        m_x[0] = 10.0 * rand() / RAND_MAX;
        m_x[1] = 10.0 * rand() / RAND_MAX;
        m_x[2] = 10.0 * rand() / RAND_MAX;
        m_t = 0.0;
    }

    inline void do_step( const double dt )
    {
        gsl_odeiv_step_apply ( m_s , m_t , dt , m_x , m_x_err , 0 , 0 , &m_sys );
        //m_t += dt;
    }

    double state( const size_t i ) const
    { return m_x[i]; }

    ~gsl_wrapper()
    {
        gsl_odeiv_step_free( m_s );
    }

private:
    double m_x[dim];
    double m_x_err[dim];
    double m_t;
    gsl_odeiv_step *m_s;
    gsl_odeiv_system m_sys;
};



int main()
{
    gsl_wrapper stepper;

    run( stepper , 20000000 / 3 , 1E-10 * 3);
}
