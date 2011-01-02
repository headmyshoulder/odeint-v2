#include <iostream>
#include <algorithm>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/timer.hpp>

#include <boost/numeric/odeint.hpp>

#define tab "\t" 

using namespace std;
using namespace boost::accumulators;
using namespace boost::numeric::odeint;

typedef accumulator_set<
    double , stats< tag::mean , tag::variance >
    > accumulator_type;

typedef boost::timer timer_type;

ostream& operator<<( ostream& out , accumulator_type &acc )
{
    out << boost::accumulators::mean( acc ) << tab;
//    out << boost::accumulators::variance( acc ) << tab;
    return out;
}

const size_t n = 3;
const double dt = 0.01;

typedef vector< double > state_type;
typedef stepper_rk4< state_type > stepper_type;

void lorenz( state_type &x , state_type &dxdt , double t )
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
    const size_t num_of_steps = 1024 * 1024 * 16;
    state_type x1( 3 , 0.0 ) , x2( 3 , 0.0 );
    stepper_type stepper;
    stepper.adjust_size( x1 );

    x1[0] = x2[0] = 1.0;

    accumulator_type acc1 , acc2;
    timer_type timer;

    size_t count = 0;
    clog.precision(4);
    while( true )
    {
        timer.restart();
	double t1 = 0.0;
	for( size_t i=0 ; i<num_of_steps ; ++i,t1+=dt )
	    stepper.do_step( lorenz , x1 , t1 , dt );
        acc1( timer.elapsed() );

        timer.restart();
	for( size_t i=0 ; i<num_of_steps; ++i )
	    stepper.do_step( lorenz , x2 , 0.0 , dt );
        acc2( timer.elapsed() );

	if( x1[0] != x2[0] ) cerr << "Error!" << endl;

	++count;

        clog << count << tab;
        clog << acc1 << tab;
        clog << acc2 << tab;

        clog << endl;
    }

    return 0;
}

