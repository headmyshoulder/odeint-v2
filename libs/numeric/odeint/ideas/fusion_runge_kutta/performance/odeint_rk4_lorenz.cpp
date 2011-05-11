/*
 * odeint_rk4_lorenz.cpp
 *
 *  Created on: May 11, 2011
 *      Author: mario
 */

#include <boost/array.hpp>

#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>
#include <boost/numeric/odeint/algebra/array_algebra.hpp>

#include "rk_performance_test_case.hpp"

typedef boost::array< double , 3 > state_type;
typedef boost::numeric::odeint::explicit_rk4< state_type , double , state_type , double ,
                                              boost::numeric::odeint::array_algebra > rk4_odeint_type;

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

inline void lorenz_func( const state_type &x , state_type &dxdt , const double t )
{
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

class odeint_wrapper
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
        m_stepper.do_step( lorenz() , m_x , m_t , dt );
        //m_t += dt;
    }

    double state( const size_t i ) const
    { return m_x[i]; }

private:
    state_type m_x;
    double m_t;
    rk4_odeint_type m_stepper;
};



int main()
{
    odeint_wrapper stepper;

    run( stepper );
}
