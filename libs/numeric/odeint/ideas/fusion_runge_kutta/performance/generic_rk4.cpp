/*
 * generic_rk4.cpp
 *
 *  Created on: Apr 28, 2011
 *      Author: mario
 */

#include <iostream>
#include <fstream>

#include <boost/array.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/timer.hpp>

//#include "fusion_explicit_rk.hpp"
#include "../fusion_explicit_rk_new.hpp"

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


typedef boost::array< double , 3 > state_type;
typedef explicit_rk< state_type , 4 > rk4_fusion_type;


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
    const double dt = 1E-10;

    accumulator_type acc;
    timer_type timer;

    srand( 12312354 );

    while( true )
    {
        state_type x = {{ 10.0 * rand()/RAND_MAX , 10.0 * rand()/RAND_MAX , 10.0 * rand()/RAND_MAX }};
        //state_type x = {{ 10.0 , 1.0 , 5.0 }};
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
