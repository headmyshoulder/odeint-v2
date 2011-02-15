/*
 * hamiltonian_steppers.cpp
 *
 *  Created on: Feb 12, 2011
 *      Author: karsten
 */

#include <iostream>
#include <tr1/array>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include <boost/timer.hpp>

#include <boost/numeric/odeint/stepper/symplectic_euler.hpp>
#include <boost/numeric/odeint/stepper/symplectic_rkn_sb3a_mclachlan.hpp>

#include <boost/numeric/odeint/hamiltonian_stepper_rk.hpp>
#include <boost/numeric/odeint/container_traits_tr1_array.hpp>


using namespace std;
using namespace boost::accumulators;


typedef accumulator_set< double , stats< tag::mean > > accumulator_type;

ostream& operator<<( ostream& out , accumulator_type &acc )
{
    out << boost::accumulators::mean( acc ) << "\t";
    return out;
}

typedef boost::timer timer_type;






typedef std::tr1::array< double , 1 > container_type;
typedef std::pair< container_type , container_type > state_type;



struct harm_osc
{
	const double m_omega_sq;

	harm_osc( double omega_sq ) : m_omega_sq( omega_sq ) { }

	template< class Coor , class MomentumDeriv >
	void operator()( const Coor &coor , MomentumDeriv &deriv )
	{
		deriv[0] = -m_omega_sq * coor[0];
	}
};


void test_euler( void )
{
	typedef boost::numeric::odeint::symplectic_euler< container_type > stepper_type;
	stepper_type stepper;
	stepper_type stepper2( stepper );
	stepper_type stepper3;
	stepper3 = stepper;

	state_type state;
	state.first[0] = 1.0;
	state.second[0] = 0.0;

	double t = 0.0;
	const double dt = 0.1;
	const double omega_sq = 4.0;

	stepper.do_step( harm_osc( omega_sq ) , state , t , dt );
}

void test_rkn_sb3a_mclachlan( void )
{
	typedef boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan< container_type > stepper_type;
	stepper_type stepper;
	stepper_type stepper2( stepper );
	stepper_type stepper3;
	stepper3 = stepper;

	state_type state;
	state.first[0] = 1.0;
	state.second[0] = 0.0;

	double t = 0.0;
	const double dt = 0.1;
	const double omega_sq = 4.0;

	stepper.do_step( harm_osc( omega_sq ) , state , t , dt );
}

void compare_euler_rkn( void )
{
	typedef boost::numeric::odeint::symplectic_euler< container_type > stepper_type1;
	typedef boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan< container_type > stepper_type2;

	stepper_type1 stepper1;
	stepper_type2 stepper2;

	state_type state1;
	state1.first[0] = 1.0;
	state1.second[0] = 0.0;
	state_type state2 = state1;

	double t = 0.0;
	const double dt = 0.1;
	const double omega_sq = 4.0;

	const size_t ilen = 100;
	for( size_t i=0 ; i<10000 ; ++i )
	{
		double q1 = state1.first[0] , p1 = state1.second[0] , q2 = state2.first[0] , p2 = state2.second[0] ;
		double energy1 = 0.5 * p1 * p1 + 0.5 * omega_sq * q1 * q1;
		double energy2 = 0.5 * p2 * p2 + 0.5 * omega_sq * q2 * q2;
		cout << t << "\t" << q1 << "\t" << p1 << "\t" << q2 << "\t" << p2 << "\t" << energy1 << "\t" << energy2 << "\n";
		for( size_t ii=0 ; ii<ilen ; ++ii , t+= dt )
		{
			stepper1.do_step( harm_osc( omega_sq ) , state1 , t , dt );
			stepper2.do_step( harm_osc( omega_sq ) , state2 , t , dt );
		}
	}
}

void performance_compare( void )
{
	typedef boost::numeric::odeint::hamiltonian_stepper_rk_qfunc< container_type > stepper_type1;
	typedef boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan< container_type > stepper_type2;

	state_type state1;
	state1.first[0] = 1.0;
	state1.second[0] = 1.0;
	state_type state2 = state1;

	const size_t steps = 100000000;

	timer_type timer;
	accumulator_type acc1 , acc2;
	harm_osc osc( 4.0 );
	const double dt = 0.01;

	for( size_t count=0 ; count<10000 ; ++count )
	{
		timer.restart();
		stepper_type1 stepper1;
		for( size_t i=0 ; i<steps ; ++i )
			stepper1.do_step( osc , state1 , 0.0 , dt );
		acc1( timer.elapsed() );

		timer.restart();
		stepper_type2 stepper2;
		for( size_t i=0 ; i<steps ; ++i )
			stepper2.do_step( osc , state2 , 0.0 , dt );
		acc2( timer.elapsed() );



		clog << count << "\t";
		clog << acc1 << "\t" << acc2 << "\t";
		clog << state1.first[0] << "\t" << state1.second[0] << "\t";
		clog << state2.first[0] << "\t" << state2.second[0] << endl;
	}


}


int main( int argc , char **argv )
{
	test_euler();
	test_rkn_sb3a_mclachlan();
	compare_euler_rkn();
	performance_compare();


	return 0;
}
