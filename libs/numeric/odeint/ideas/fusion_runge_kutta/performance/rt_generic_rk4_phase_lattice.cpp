/*
 * rt_generic_rk4_phase_lattice.cpp
 *
 *  Created on: May 15, 2011
 *      Author: mario
 */

#include <boost/array.hpp>

#include "../rt_explicit_rk.hpp"

#include "rk_performance_test_case.hpp"

#include "phase_lattice.hpp"

const size_t N = 1024;

typedef boost::array< double , N > state_type;

typedef rt_explicit_rk< state_type > rk_stepper_type;

const size_t stage_count = 4;


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
        for( size_t i = 0 ; i<N ; ++i )
			m_x[i] = 2.0*3.1415927*rand() / RAND_MAX;
        m_t = 0.0;
    }

    inline void do_step( const double dt )
    {
        m_stepper.do_step( phase_lattice<N>() , m_x , m_t , dt );
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
	srand( 12312354 );

    rt_generic_wrapper stepper;

    run( stepper , 10000 , 1E-6 );
}
