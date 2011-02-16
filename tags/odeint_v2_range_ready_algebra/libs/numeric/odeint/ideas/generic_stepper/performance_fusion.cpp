/*
 * performance.cpp
 *
 *  Created on: Dec 1, 2010
 *      Author: mario
 */

/*
 * butcher_test.cpp
 *
 *  Created on: Nov 5, 2010
 *      Author: karsten
 */

#include <iostream>
#include <fstream>

#include <tr1/array>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/timer.hpp>

#include "runge_kutta_stepper_fast.hpp"

#define tab "\t"

using namespace std;
using namespace boost::accumulators;

typedef accumulator_set<
    double , stats< tag::mean , tag::variance >
    > accumulator_type;

ostream& operator<<( ostream& out , accumulator_type &acc )
{
    out << boost::accumulators::mean( acc ) << tab;
//    out << boost::accumulators::variance( acc ) << tab;
    return out;
}

typedef boost::timer timer_type;


typedef std::tr1::array< double , 3 > state_type;
typedef runge_kutta_stepper< state_type , 4 > rk4_fusion_type;


void lorenz( const state_type &x , state_type &dxdt , double t )
{
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}




int main( int argc , char **argv )
{
    typedef rk4_fusion_type::coef_a_type coef_a_type;
    typedef rk4_fusion_type::coef_b_type coef_b_type;
    typedef rk4_fusion_type::coef_c_type coef_c_type;

    const boost::array< double , 1 > a1 = {{ 0.5 }};
    const boost::array< double , 2 > a2 = {{ 0.0 , 0.5 }};
    const boost::array< double , 3 > a3 = {{ 0.0 , 0.0 , 1.0 }};

    const coef_a_type a = fusion::make_vector( a1 , a2 , a3 );
    const coef_b_type b = {{ 1.0/6 , 1.0/3 , 1.0/3 , 1.0/6 }};
    const coef_c_type c = {{ 0.0 , 0.5 , 0.5 , 1.0 }};

    rk4_fusion_type rk4_fusion( a , b , c );

    const size_t num_of_steps = 20000000;
    const size_t dt = 0.01;

    accumulator_type acc;
    timer_type timer;

    srand48( 12312354 );

    while( true )
    {
        state_type x = {{ 10.0 * drand48() , 10.0 * drand48() , 10.0 * drand48() }};
        double t = 0.0;

        timer.restart();
        for( size_t i=0 ; i<num_of_steps ; ++i, t+=dt )
            rk4_fusion.do_step( lorenz , x , t , dt );
        acc( timer.elapsed() );

        clog.precision( 3 );
        clog.width( 5 );
        clog << acc << " " << x[0] << endl;
    }



    return 0;
}

/*
 * Compile with :
 * g++ -O3 -I$BOOST_ROOT -I$HOME/boost/chrono -I$ODEINT_ROOT butcher_performance.cpp
 */
