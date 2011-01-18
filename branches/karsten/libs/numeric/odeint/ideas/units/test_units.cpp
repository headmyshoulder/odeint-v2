/*
 * test.cpp
 *
 *  Created on: Jan 18, 2011
 *      Author: karsten
 */

#include <iostream>
#include <tr1/array>

#include <boost/units/systems/si/energy.hpp>
#include <boost/units/systems/si/force.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/electric_potential.hpp>
#include <boost/units/systems/si/current.hpp>
#include <boost/units/systems/si/resistance.hpp>
#include <boost/units/systems/si/io.hpp>

#include "explicit_euler_units.hpp"

namespace units = boost::units;
namespace si = boost::units::si;

using std::cout;
using std::endl;

typedef units::derived_dimension< units::length_base_dimension , 1 , units::time_base_dimension , -1 > vel_type;

typedef std::tr1::array< units::quantity< si::length > , 2 > state_type;
typedef std::tr1::array< units::quantity< vel_type > , 2 > deriv_type;
typedef units::quantity< si::time > time_type;

typedef boost::numeric::odeint::explicit_euler_units< deriv_type , double , time_type > stepper_type;

struct lorenz
{
	template< class State , class Deriv , class Time >
	void operator()( const State &x , Deriv &dxdt , const Time &t )
	{
		dxdt[0] = ( x[0] / Time( 1.0 * si::second ) );
		dxdt[0] = ( x[1] / Time( 1.0 * si::second ) );
	}
};



int main( int argc , char **argv )
{
	/*
	 * zwei faelle, erstens units, zweitens ranges
	 */
	state_type x = {{ 1.0 * si::meter , 1.0 * si::meter }};
	time_type t( 0.0 * si::seconds ) , dt( 0.0 * si::seconds );
	stepper_type stepper;
	stepper.do_step( lorenz() , x , t , dt );

//	units::quantity< si::length > x( 2.0 * si::meter );
//	units::quantity< si::time > t( 1.0 * si::second );

//	cout << x * t << endl;

	return 0;
}
