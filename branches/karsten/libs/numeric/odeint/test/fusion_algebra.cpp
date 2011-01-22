/*
 * check_operations.cpp
 *
 *  Created on: Jan 20, 2011
 *      Author: karsten
 */

#define BOOST_TEST_MODULE odeint_fusion_algebra
#include <cmath>
#include <complex>
#include <utility>
#include <functional>
#include <tr1/array>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <boost/range.hpp>

#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/velocity.hpp>
#include <boost/units/systems/si/acceleration.hpp>
#include <boost/units/systems/si/io.hpp>

#include <boost/fusion/container.hpp>

#include <boost/numeric/odeint/algebra/standard_operations.hpp>
#include <boost/numeric/odeint/algebra/fusion_algebra.hpp>

namespace units = boost::units;
namespace si = boost::units::si;
namespace fusion = boost::fusion;

using boost::numeric::odeint::standard_operations;
using boost::numeric::odeint::fusion_algebra;

typedef double value_type;
typedef units::quantity< si::time , value_type > time_type;
typedef units::quantity< si::length , value_type > length_type;
typedef units::quantity< si::velocity , value_type > velocity_type;
typedef units::quantity< si::acceleration , value_type > acceleration_type;
typedef fusion::vector< length_type , velocity_type > state_type;
typedef fusion::vector< velocity_type , acceleration_type > deriv_type;

const time_type dt = 0.1 * si::second;


struct fusion_fixture
{
	fusion_fixture( void )
	: res( 0.0 * si::meter , 0.0 * si::meter_per_second ) ,
	  x( 1.0 * si::meter , 1.0 * si::meter_per_second ) ,
	  k1( 1.0 * si::meter_per_second , 1.0 * si::meter_per_second_squared ) ,
	  k2( 2.0 * si::meter_per_second , 2.0 * si::meter_per_second_squared ) ,
	  k3( 3.0 * si::meter_per_second , 3.0 * si::meter_per_second_squared ) ,
	  k4( 4.0 * si::meter_per_second , 4.0 * si::meter_per_second_squared ) ,
	  k5( 5.0 * si::meter_per_second , 5.0 * si::meter_per_second_squared ) ,
	  k6( 6.0 * si::meter_per_second , 6.0 * si::meter_per_second_squared ) { }

	~fusion_fixture( void )
	{
		BOOST_CHECK_CLOSE( fusion::at_c< 0 >( x ).value() , 1.0 , 1.0e-10 );
		BOOST_CHECK_CLOSE( fusion::at_c< 1 >( x ).value() , 1.0 , 1.0e-10 );
		BOOST_CHECK_CLOSE( fusion::at_c< 0 >( k1 ).value() , 1.0 , 1.0e-10 );
		BOOST_CHECK_CLOSE( fusion::at_c< 1 >( k1 ).value() , 1.0 , 1.0e-10 );
		BOOST_CHECK_CLOSE( fusion::at_c< 0 >( k2 ).value() , 2.0 , 1.0e-10 );
		BOOST_CHECK_CLOSE( fusion::at_c< 1 >( k2 ).value() , 2.0 , 1.0e-10 );
		BOOST_CHECK_CLOSE( fusion::at_c< 0 >( k3 ).value() , 3.0 , 1.0e-10 );
		BOOST_CHECK_CLOSE( fusion::at_c< 1 >( k3 ).value() , 3.0 , 1.0e-10 );
		BOOST_CHECK_CLOSE( fusion::at_c< 0 >( k4 ).value() , 4.0 , 1.0e-10 );
		BOOST_CHECK_CLOSE( fusion::at_c< 1 >( k4 ).value() , 4.0 , 1.0e-10 );
		BOOST_CHECK_CLOSE( fusion::at_c< 0 >( k5 ).value() , 5.0 , 1.0e-10 );
		BOOST_CHECK_CLOSE( fusion::at_c< 1 >( k5 ).value() , 5.0 , 1.0e-10 );
		BOOST_CHECK_CLOSE( fusion::at_c< 0 >( k6 ).value() , 6.0 , 1.0e-10 );
		BOOST_CHECK_CLOSE( fusion::at_c< 1 >( k6 ).value() , 6.0 , 1.0e-10 );
	}

	state_type res , x;
	deriv_type k1 , k2 , k3 , k4 , k5 , k6 ;
};



BOOST_AUTO_TEST_SUITE( fusion_algebra_test )

BOOST_AUTO_TEST_CASE( for_each2 )
{
	fusion_fixture f;
	fusion_algebra::for_each2()( f.res , f.k1 ,
			standard_operations::scale_sum1< time_type >( dt ) );
	BOOST_CHECK_CLOSE( fusion::at_c< 0 >( f.res ).value() , 0.1 , 1.0e-10 );
	BOOST_CHECK_CLOSE( fusion::at_c< 1 >( f.res ).value() , 0.1 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( for_each3 )
{
	fusion_fixture f;
	fusion_algebra::for_each3()( f.res , f.x , f.k1 ,
			standard_operations::scale_sum2< value_type , time_type >( 1.0 , dt ) );
	BOOST_CHECK_CLOSE( fusion::at_c< 0 >( f.res ).value() , 1.1 , 1.0e-10 );
	BOOST_CHECK_CLOSE( fusion::at_c< 1 >( f.res ).value() , 1.1 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( for_each4 )
{
	fusion_fixture f;
	fusion_algebra::for_each4()( f.res , f.x , f.k1 , f.k2 ,
			standard_operations::scale_sum3< value_type , time_type , time_type >( 1.0 , dt , dt ) );
	BOOST_CHECK_CLOSE( fusion::at_c< 0 >( f.res ).value() , 1.3 , 1.0e-10 );
	BOOST_CHECK_CLOSE( fusion::at_c< 1 >( f.res ).value() , 1.3 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( for_each5 )
{
	fusion_fixture f;
	fusion_algebra::for_each5()( f.res , f.x , f.k1 , f.k2 , f.k3 ,
			standard_operations::scale_sum4< value_type , time_type , time_type , time_type >( 1.0 , dt , dt , dt ) );
	BOOST_CHECK_CLOSE( fusion::at_c< 0 >( f.res ).value() , 1.6 , 1.0e-10 );
	BOOST_CHECK_CLOSE( fusion::at_c< 1 >( f.res ).value() , 1.6 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( for_each6 )
{
	fusion_fixture f;
	fusion_algebra::for_each6()( f.res , f.x , f.k1 , f.k2 , f.k3 , f.k4 ,
			standard_operations::scale_sum5< value_type , time_type , time_type , time_type , time_type >( 1.0 , dt , dt , dt , dt ) );
	BOOST_CHECK_CLOSE( fusion::at_c< 0 >( f.res ).value() , 2.0 , 1.0e-10 );
	BOOST_CHECK_CLOSE( fusion::at_c< 1 >( f.res ).value() , 2.0 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( for_each7 )
{
	fusion_fixture f;
	fusion_algebra::for_each7()( f.res , f.x , f.k1 , f.k2 , f.k3 , f.k4 , f.k5 ,
			standard_operations::scale_sum6< value_type , time_type , time_type , time_type , time_type , time_type >( 1.0 , dt , dt , dt , dt , dt ) );
	BOOST_CHECK_CLOSE( fusion::at_c< 0 >( f.res ).value() , 2.5 , 1.0e-10 );
	BOOST_CHECK_CLOSE( fusion::at_c< 1 >( f.res ).value() , 2.5 , 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( for_each8 )
{
	fusion_fixture f;
	fusion_algebra::for_each8()( f.res , f.x , f.k1 , f.k2 , f.k3 , f.k4 , f.k5 , f.k6 ,
			standard_operations::scale_sum7< value_type , time_type , time_type , time_type , time_type , time_type , time_type >( 1.0 , dt , dt , dt , dt , dt ,dt ) );
	BOOST_CHECK_CLOSE( fusion::at_c< 0 >( f.res ).value() , 3.1 , 1.0e-10 );
	BOOST_CHECK_CLOSE( fusion::at_c< 1 >( f.res ).value() , 3.1 , 1.0e-10 );
}


BOOST_AUTO_TEST_CASE( reduce )
{
	double sum = fusion_algebra::reduce()( fusion::make_vector( 1.0 , 2.0 , 3.0 ) , std::plus< double >() , 0.0 );
	BOOST_CHECK_CLOSE( sum , 6.0 , 1.0e-10 );
}


BOOST_AUTO_TEST_SUITE_END()
