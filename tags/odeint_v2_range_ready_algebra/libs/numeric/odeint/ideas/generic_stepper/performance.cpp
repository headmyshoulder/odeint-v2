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

#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/timer.hpp>

#include "predefined_steppers.hpp"
#include "runge_kutta_stepper.hpp"

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
typedef boost::numeric::odeint::explicit_rk4< state_type > rk4_odeint_type;
typedef mpl_rk4_stepper< state_type > rk4_mpl_type;
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
    rk4_odeint_type rk4_odeint;
    rk4_mpl_type rk4_mpl;

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

    accumulator_type acc1 , acc2 , acc3;
    timer_type timer;

    srand48( 12312354 );

    while( true )
    {
        state_type x1 = {{ 10.0 * drand48() , 10.0 * drand48() , 10.0 * drand48() }};
        state_type x2 = x1;
        state_type x3 = x1;
        double t1 = 0.0 ;
        double t2 = 0.0;
        double t3 = 0.0;

        timer.restart();
        for( size_t i=0 ; i<num_of_steps ; ++i,t1+=dt )
            rk4_odeint.do_step( lorenz , x1 , t1 , dt );
        acc1( timer.elapsed() );

        timer.restart();
        for( size_t i=0 ; i<num_of_steps ; ++i,t2+=dt )
            rk4_mpl.do_step( lorenz , x2 , t2 , dt );
        acc2( timer.elapsed() );

        timer.restart();
        for( size_t i=0 ; i<num_of_steps ; ++i,t3+=dt )
            rk4_fusion.do_step( lorenz , x3 , t3 , dt );
        acc3( timer.elapsed() );

        clog.precision( 3 );
        clog.width( 5 );
        clog << acc1 << " " << acc2 << " " << acc3 << " " << x1[0] << " " << x2[0] << " " << x3[0] << endl;
    }



    return 0;
}

/*
 * Compile with :
 * g++ -O3 -I$BOOST_ROOT -I$HOME/boost/chrono -I$ODEINT_ROOT butcher_performance.cpp
 */
