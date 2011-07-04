/*
 * generic_rk54ck_lorenz.cpp
 *
 *  Created on: May 12, 2011
 *      Author: mario
 */

#include <boost/array.hpp>

#include "../fusion_explicit_error_rk.hpp"

#include "rk_performance_test_case.hpp"

#include "lorenz.hpp"

typedef boost::array< double , 3 > state_type;
typedef explicit_error_rk< state_type , 6 > rk54ck_fusion_type;


typedef rk54ck_fusion_type::coef_a_type coef_a_type;
typedef rk54ck_fusion_type::coef_b_type coef_b_type;
typedef rk54ck_fusion_type::coef_c_type coef_c_type;

const boost::array< double , 1 > a1 = {{ 0.2 }};
const boost::array< double , 2 > a2 = {{ 3.0/40.0 , 9.0/40 }};
const boost::array< double , 3 > a3 = {{ 0.3 , -0.9 , 1.2 }};
const boost::array< double , 4 > a4 = {{ -11.0/54.0 , 2.5 , -70.0/27.0 , 35.0/27.0 }};
const boost::array< double , 5 > a5 = {{ 1631.0/55296.0 , 175.0/512.0 , 575.0/13824.0 , 44275.0/110592.0 , 253.0/4096.0 }};

const coef_a_type a = fusion::make_vector( a1 , a2 , a3 , a4 , a5 );
const coef_b_type b = {{ 37.0/378.0 , 0.0 , 250.0/621.0 , 125.0/594.0 , 0.0 , 512.0/1771.0 }};
const coef_b_type b2 = {{ b[0]-2825.0/27648.0 , b[1]-0.0 , b[2]-18575.0/48384.0 , b[3]-13525.0/55296.0 , b[4]-277.0/14336.0 , b[5]-0.25 }};

const coef_c_type c = {{ 0.0 , 0.2 , 0.3 , 0.6 , 1.0 , 7.0/8 }};

class fusion_wrapper
{
public:

    fusion_wrapper() : m_stepper( a , b , b2 , c )
    { }

    void reset_init_cond()
    {
        m_x[0] = 10.0 * rand() / RAND_MAX;
        m_x[1] = 10.0 * rand() / RAND_MAX;
        m_x[2] = 10.0 * rand() / RAND_MAX;
        m_t = 0.0;
    }

    inline void do_step( const double dt )
    {
        m_stepper.do_step( lorenz() , m_x , m_t , dt , m_x_err );
        //m_t += dt;
    }

    double state( const size_t i ) const
    { return m_x[i]; }

private:
    state_type m_x;
    state_type m_x_err;
    double m_t;
    rk54ck_fusion_type m_stepper;
};



int main()
{
    srand( 12312354 );

    fusion_wrapper stepper;

    run( stepper );
}
