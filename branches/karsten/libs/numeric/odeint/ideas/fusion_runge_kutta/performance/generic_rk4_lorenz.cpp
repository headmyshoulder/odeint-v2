#include <boost/array.hpp>

#include "../fusion_explicit_rk_new.hpp"

#include "rk_performance_test_case.hpp"

typedef boost::array< double , 3 > state_type;
typedef explicit_rk< state_type , 4 > rk4_fusion_type;

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
        m_x[0] = 10.0 * rand() / RAND_MAX;
        m_x[1] = 10.0 * rand() / RAND_MAX;
        m_x[2] = 10.0 * rand() / RAND_MAX;
        m_t = 0.0;
    }

    inline void do_step( const double dt )
    {
        m_stepper.do_step( lorenz() , m_x , m_t , dt );
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
    fusion_wrapper stepper;

    run( stepper );
}
