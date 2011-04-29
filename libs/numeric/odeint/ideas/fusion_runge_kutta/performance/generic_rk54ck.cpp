/*
 * generic_rk78.cpp
 *
 *  Created on: Apr 29, 2011
 *      Author: mario
 */


#include <iostream>
#include <fstream>

#include <boost/array.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/timer.hpp>

#include "../fusion_explicit_error_rk.hpp"

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
typedef explicit_error_rk< state_type , 6 > rk54ck_fusion_type;


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
    typedef rk54ck_fusion_type::coef_a_type coef_a_type;
    typedef rk54ck_fusion_type::coef_b_type coef_b_type;
    typedef rk54ck_fusion_type::coef_c_type coef_c_type;

    const boost::array< double , 1 > a1 = {{ 0.2 }};
    const boost::array< double , 2 > a2 = {{ 3.0/40 , 9.0/40 }};
    const boost::array< double , 3 > a3 = {{ 3.0/10 , -9.0/10 , 6.0/5 }};
    const boost::array< double , 4 > a4 = {{ -11.0/54 , 5.0/2 , -70.0/27 , 35.0/27 }};
    const boost::array< double , 5 > a5 = {{ 1631.0/55296 , 175.0/512 , 575.0/13824 , 44275.0/110592 , 253.0/4096 }};

    const coef_a_type a = fusion::make_vector( a1 , a2 , a3 , a4 , a5 );
    const coef_b_type b = {{ 37.0/378 , 0.0 , 250.0/621 , 125.0/594 , 0.0 , 512.0/1771 }};
    const coef_b_type b2 = {{ b[0]-2825.0/27648 , b[1]-0.0 , b[2]-18575.0/48384 , b[3]-13525.0/55296 , b[4]-277.0/14336 , b[5]-1.0/4 }};
    //const coef_b_type b = {{ 2825.0/27648 , 0.0 , 18575.0/48384 , 13525.0/55296 , 277.0/14336 , 1.0/4 }};
    //const coef_b_type b2 = {{ b[0]-37.0/378 , b[1]-0.0 , b[2]-250.0/621 , b[3]-125.0/594 , b[4]-0.0 , b[5]-512.0/1771 }};
    const coef_c_type c = {{ 0.0 , 0.2 , 0.3 , 0.6 , 1.0 , 7.0/8 }};

    rk54ck_fusion_type rk54ck( a , b , b2 , c );

    const size_t num_of_steps = 20000000;
    const double dt = 1E-4;

    accumulator_type acc;
    timer_type timer;

    srand( 12312354 );

    while( true )
    {
        state_type x = {{ 10.0 * rand()/RAND_MAX , 10.0 * rand()/RAND_MAX , 10.0 * rand()/RAND_MAX }};
        state_type x_err;
        double t = 0.0;

        timer.restart();
        for( size_t i=0 ; i<num_of_steps ; ++i, t+=dt )
            rk54ck.do_step( lorenz , x , t , dt , x_err );
        acc( timer.elapsed() );

        clog.precision( 15 );
        clog.width( 20 );
        clog << acc << " " << x[0] << tab << " " << x_err[0] << endl;
    }

    return 0;
}
