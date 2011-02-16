/*
 * prepare_stepper_testing.hpp
 *
 *  Created on: Jan 19, 2011
 *      Author: karsten
 */

#ifndef PREPARE_STEPPER_TESTING_HPP_
#define PREPARE_STEPPER_TESTING_HPP_

#include <tr1/array>
#include <vector>

#include <boost/mpl/vector.hpp>

namespace mpl = boost::mpl;


struct constant_system_standard
{
	template< class State , class Deriv , class Time >
	void operator()( const State &x , Deriv &dxdt , const Time &t ) const
	{
		dxdt[0] = 1.0;
	}
};

struct constant_system_vector_space
{
	template< class State , class Deriv , class Time >
	void operator()( const State &x , class Deriv &dxdt , const Time &t  ) const
	{
		dxdt.m_x = 1.0;
	}
};

struct constant_system_fusion
{
	template< class State , class Deriv , class Time >
	void operator()( const State &x , class Deriv &dxdt , const Time &t ) const
	{
		fusion::at_c< 0 >( dxdt ) = fusion::at_c< 0 >( x ) / Time( 1.0 );
	}
};



typedef mpl::vector
<
	mpl::vector< float , std::tr1::array< float , 1 > , std::tr1::array< float , 1 > , constant_system_standard , stepper_type > ,
    mpl::vector< float , std::tr1::array< float , 1 > , std::vector< float > , constant_system_standard , stepper_type >
> types_and_systems_matrix;


struct check_stepper
{
	template< class Vector >
	void operator()( Vector )
	{
		typedef typename mpl::at_c< Vector , 0 >::type value_type;
		typedef typename mpl::at_c< Vector , 1 >::type deriv_type;
		typedef typename mpl::at_c< Vector , 2 >::type state_type;
		typedef typename mpl::at_c< Vector , 3 >::type system_type;
	}
};


#endif /* PREPARE_STEPPER_TESTING_HPP_ */
