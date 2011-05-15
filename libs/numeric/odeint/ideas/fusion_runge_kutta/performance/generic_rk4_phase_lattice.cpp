/*
 * generic_rk4_phase_lattice.cpp
 *
 *  Created on: May 15, 2011
 *      Author: mario
 */

#include <boost/array.hpp>

#include "../fusion_explicit_rk_new.hpp"

#include "rk_performance_test_case.hpp"

#include "phase_lattice.hpp"

const size_t N = 1024;

typedef boost::array< double , N > state_type;
typedef explicit_rk< state_type , 4 > rk4_fusion_type;

typedef rk4_fusion_type::coef_a_type coef_a_type;
typedef rk4_fusion_type::coef_b_type coef_b_type;
typedef rk4_fusion_type::coef_c_type coef_c_type;

const boost::array< double , 1 > a1 = {{ 0.5 }};
const boost::array< double , 2 > a2 = {{ 0.0 , 0.5 }};
const boost::array< double , 3 > a3 = {{ 0.0 , 0.0 , 1.0 }};

const coef_a_type a = fusion::make_vector( a1 , a2 , a3 );
const coef_b_type b = {{ 1.0/6 , 1.0/3 , 1.0/3 , 1.0/6 }};
const coef_c_type c = {{ 0.0 , 0.5 , 0.5 , 1.0 }};

class fusion_wrapper
{

public:

    fusion_wrapper() : m_stepper( a , b , c )
    { }

    void reset_init_cond()
    {
        for( size_t i = 0 ; i<N ; ++i )
			m_x[i] = 2.0*3.1415927*rand() / RAND_MAX;
		m_t = 0.0;
    }

    inline void do_step( const double dt )
    {
        m_stepper.do_step( phase_lattice<N>(), m_x , m_t , dt );
    }

    double state( const size_t i ) const
    { return m_x[i]; }

private:
    state_type m_x;
    double m_t;
    rk4_fusion_type m_stepper;
};



int main()
{
    srand( 12312354 );

    fusion_wrapper stepper;

    run( stepper , 1000 , 1E-2 );
}
