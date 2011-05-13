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

#include "lorenz.hpp"

typedef boost::array< double , 3 > state_type;
typedef boost::numeric::odeint::explicit_rk4< state_type , double , state_type , double ,
                                              boost::numeric::odeint::array_algebra > rk4_odeint_type;


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
    srand( 12312354 );

    odeint_wrapper stepper;

    run( stepper );
}
