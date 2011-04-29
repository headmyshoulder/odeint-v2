/*
 * rt_generic_rk4.cpp
 *
 *  Created on: Apr 29, 2011
 *      Author: mario
 */

#include <iostream>
#include <fstream>
#include <vector>

#include <boost/array.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/timer.hpp>

#include "rt_algebra.hpp"

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


template< class StateType >
class rt_explicit_rk
{
public:
	typedef StateType state_type;
	typedef vector< vector< double > > coeff_a_type;
	typedef vector< double > coeff_b_type;
	typedef vector< double > coeff_c_type;

	rt_explicit_rk( size_t stage_count , coeff_a_type &a , coeff_b_type &b , coeff_c_type &c )
		: m_s( stage_count ) , m_a( a ) , m_b( b ) , m_c( c )
	{ 
		m_F = new state_type[ m_s ];
	}

	~rt_explicit_rk() { delete[] m_F; }

	template< class System >
	void do_step( System &sys , state_type &x , const double t , const double dt )
	{
		// first stage separately
		sys( x , m_F[0] , t + m_c[0]*t );
		if( m_s == 1 )
			rt_algebra::foreach( x , x , m_b , m_F , dt );
		else
			rt_algebra::foreach( m_x_tmp , x , m_a[0] , m_F , dt );

		for( size_t stage = 2 ; stage <= m_s ; ++stage )
		{
			sys( m_x_tmp , m_F[stage-1] , t + m_c[stage-1]*dt );
			if( stage == m_s )
				rt_algebra::foreach( x , x , m_b , m_F , dt );
			else
				rt_algebra::foreach( m_x_tmp , x , m_a[stage-1] , m_F , dt );
		}
	}



private:
	size_t m_s;
	vector< vector< double > > m_a;
	vector< double > m_b;
	vector< double > m_c;

	state_type m_x_tmp;
	state_type *m_F;
};

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

typedef rt_explicit_rk< state_type > rk_stepper_type;

const size_t stage_count = 4;

int main( int argc , char **argv )
{
	rk_stepper_type::coeff_a_type a( stage_count-1 );
	a[0].resize(1); a[0][0] = 0.5;
	a[1].resize(2); a[1][0] = 0.0; a[1][1] = 0.5;
	a[2].resize(3); a[2][0] = 0.0; a[2][1] = 0.0; a[2][2] = 1.0;

	rk_stepper_type::coeff_b_type b( stage_count );
	b[0] = 1.0/6; b[1] = 1.0/3; b[2] = 1.0/3; b[3] = 1.0/6;

	rk_stepper_type::coeff_c_type c( stage_count );
	c[0] = 0.0; c[1] = 0.5; c[2] = 0.5; c[3] = 1.0;

    rk_stepper_type rk4_rt( stage_count , a , b , c );

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
            rk4_rt.do_step( lorenz , x , t , dt );
        acc( timer.elapsed() );

        clog.precision( 3 );
        clog.width( 5 );
        clog << acc << " " << x[0] << endl;
    }

    return 0;
}