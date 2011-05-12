/*
 * rt_generic_rk4_lorenz.cpp
 *
 *  Created on: May 11, 2011
 *      Author: mario
 */

#include <boost/array.hpp>

#include "../rt_explicit_rk.hpp"

#include "rk_performance_test_case.hpp"

typedef boost::array< double , 3 > state_type;

typedef rt_explicit_rk< state_type > rk_stepper_type;

const size_t stage_count = 4;

inline void lorenz( const state_type &x , state_type &dxdt , const double t )
{
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}



class rt_generic_wrapper
{
public:

    rt_generic_wrapper() : m_stepper( stage_count )
    {
        rk_stepper_type::coeff_a_type a( stage_count-1 );
        a[0].resize(1); a[0][0] = 0.5;
        a[1].resize(2); a[1][0] = 0.0; a[1][1] = 0.5;
        a[2].resize(3); a[2][0] = 0.0; a[2][1] = 0.0; a[2][2] = 1.0;

        rk_stepper_type::coeff_b_type b( stage_count );
        b[0] = 1.0/6; b[1] = 1.0/3; b[2] = 1.0/3; b[3] = 1.0/6;

        rk_stepper_type::coeff_c_type c( stage_count );
        c[0] = 0.0; c[1] = 0.5; c[2] = 0.5; c[3] = 1.0;

        m_stepper.set_params( a , b , c );
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
        m_stepper.do_step( lorenz , m_x , m_t , dt );
        //m_t += dt;
    }

    double state( const size_t i ) const
    { return m_x[i]; }

private:
    state_type m_x;
    double m_t;
    rk_stepper_type m_stepper;
};



int main()
{
    rt_generic_wrapper stepper;

    run( stepper );
}
