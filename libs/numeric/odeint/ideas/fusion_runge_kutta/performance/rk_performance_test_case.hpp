/*
 * rk_performance_test_case.hpp
 *
 *  Created on: May 11, 2011
 *      Author: mario
 */

#include <iostream>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/timer.hpp>

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


template< class Stepper >
void run( Stepper &stepper , const size_t num_of_steps = 20000000 , const double dt = 1E-10 )
{
    const size_t loops = 20;

    accumulator_type acc;
    timer_type timer;

    srand( 12312354 );

    for( size_t n=0 ; n<loops ; ++n )
    {
        stepper.reset_init_cond( );

        timer.restart();
        for( size_t i = 0 ; i < num_of_steps ; ++i )
            stepper.do_step( dt );
        acc(timer.elapsed());

        clog.precision(3);
        clog.width(5);
        clog << acc << " " << stepper.state(0) << endl;
    }
    cout << acc << endl;
}
