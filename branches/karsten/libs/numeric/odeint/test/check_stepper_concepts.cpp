/* Boost stepper_euler.cpp test file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file tests the use of the all different steppers with several state types:
 std::vector< double >
 gsl_vector
 vector_space_1d< double >  (see vector_space_1d.hpp)
 std::tr1::array< double , 1 >
  
 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#define BOOST_TEST_MODULE test_stepper_concepts

#include <vector>
#include <cmath>
#include <iostream>
#include <tr1/array>

#include <boost/test/unit_test.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/mpl/end.hpp>

#include <boost/bind.hpp>
#include <boost/utility.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/algebra/external/gsl_vector_adaptor.hpp>

#include "vector_space_1d.hpp"


using std::vector;

using namespace boost::unit_test;
using namespace boost::numeric::odeint;
namespace mpl = boost::mpl;


typedef std::vector< double > vector_type;
typedef gsl_vector gsl_vector_type;
typedef vector_space_1d< double > vector_space_type;
typedef std::tr1::array< double , 1 > array_type;

template< class State > struct algebra_dispatcher { typedef standard_algebra< State > type; };
template<> struct algebra_dispatcher< vector_space_type > { typedef vector_space_algebra type; };


void constant_system1( const vector_type &x , vector_type &dxdt , double t ) { dxdt[0] = 1.0; }
void constant_system2( const gsl_vector_type &x , gsl_vector_type &dxdt , double t ) { gsl_vector_set( &dxdt , 0 , 1.0 ); }
void constant_system3( const vector_space_type &x , vector_space_type &dxdt , double t ) { dxdt.m_x = 1.0; }
void constant_system4( const array_type &x , array_type &dxdt , double t ) { dxdt[0] = 1.0; }




const double eps = 1.0e-14;

template< class Stepper , class System >
void check_stepper_concept( Stepper &stepper , System system , typename Stepper::state_type &x )
{
    typedef Stepper stepper_type;
    typedef typename stepper_type::state_type container_type;
    typedef typename stepper_type::order_type order_type;
    typedef typename stepper_type::time_type time_type;

    stepper.do_step( system , x , 0.0 , 0.1 );
}

template< class Stepper , class System >
void check_error_stepper_concept( Stepper &stepper , System system , typename Stepper::state_type &x , typename Stepper::state_type &xerr )
{
    typedef Stepper stepper_type;
    typedef typename stepper_type::state_type container_type;
    typedef typename stepper_type::order_type order_type;
    typedef typename stepper_type::time_type time_type;

    stepper.do_step( system , x , 0.0 , 0.1 , xerr );
}

template< class Stepper , class System >
void check_controlled_stepper_concept( Stepper &stepper , System system , typename Stepper::state_type &x )
{
	typedef Stepper stepper_type;
    typedef typename stepper_type::state_type container_type;
    typedef typename stepper_type::order_type order_type;
    typedef typename stepper_type::time_type time_type;

    time_type t = 0.0 , dt = 0.1;
    stepper.try_step( system , x , t , dt );
}





template< class Stepper , class State > struct perform_stepper_test;

template< class Stepper >
struct perform_stepper_test< Stepper , vector_type >
{
	void operator()( void )
	{
		vector_type x( 1 , 0.0 );
		Stepper stepper;
		check_stepper_concept( stepper , constant_system1 , x );
		BOOST_CHECK_SMALL( fabs( x[0] - 0.1 ) , eps );
	}
};

template< class Stepper >
struct perform_stepper_test< Stepper , gsl_vector_type >
{
	void operator()( void ) const
	{
		gsl_vector_type *x = gsl_vector_alloc( 1 );
		gsl_vector_set( x , 0 , 0.0 );
		Stepper stepper;
		check_stepper_concept( stepper , constant_system2 , *x );
		BOOST_CHECK_SMALL( fabs( gsl_vector_get( x , 0 ) - 0.1 ) , eps );
		gsl_vector_free( x );
	}
};

template< class Stepper >
struct perform_stepper_test< Stepper , vector_space_type >
{
	void operator()( void ) const
	{
		vector_space_type x;
		x.m_x = 0.0;
		Stepper stepper;
		check_stepper_concept( stepper , constant_system3 , x );
		BOOST_CHECK_SMALL( fabs( x.m_x - 0.1 ) , eps );
	}
};

template< class Stepper >
struct perform_stepper_test< Stepper , array_type >
{
	void operator()( void )
	{
		array_type x;
		x[0] = 0.0;
		Stepper stepper;
		check_stepper_concept( stepper , constant_system4 , x );
		BOOST_CHECK_SMALL( fabs( x[0] - 0.1 ) , eps );
	}
};




template< class State > class stepper_methods : public mpl::vector<
	explicit_euler< State , double , typename algebra_dispatcher< State >::type > ,
	explicit_rk4< State , double , typename algebra_dispatcher< State >::type > ,
	explicit_error_rk54_ck< State , double , typename algebra_dispatcher< State >::type >
> { };


typedef mpl::insert_range<	mpl::vector0<> , mpl::end< mpl::vector0<> >::type ,	stepper_methods< vector_type > >::type first_stepper_type;
typedef mpl::insert_range<	first_stepper_type , mpl::end< first_stepper_type >::type , stepper_methods< gsl_vector_type > >::type second_stepper_type;
typedef mpl::insert_range<	second_stepper_type , mpl::end< second_stepper_type >::type , stepper_methods< vector_space_type > >::type third_stepper_type;
typedef mpl::insert_range<	third_stepper_type , mpl::end< third_stepper_type >::type , stepper_methods< array_type > >::type all_stepper_methods;



BOOST_AUTO_TEST_CASE_TEMPLATE( stepper_test , Stepper, all_stepper_methods )
{
	perform_stepper_test< Stepper , typename Stepper::state_type > tester;
	tester();
}




















template< class Stepper , class State > struct perform_error_stepper_test;

template< class Stepper >
struct perform_error_stepper_test< Stepper , vector_type >
{
	void operator()( void )
	{
		vector_type x( 1 , 0.0 ) , xerr( 1 );
		Stepper stepper;
		check_error_stepper_concept( stepper , constant_system1 , x , xerr );
		BOOST_CHECK_SMALL( fabs( x[0] - 0.1 ) , eps );
	}
};

template< class Stepper >
struct perform_error_stepper_test< Stepper , gsl_vector_type >
{
	void operator()( void ) const
	{
		gsl_vector_type *x = gsl_vector_alloc( 1 ) , *xerr = gsl_vector_alloc( 1 );
		gsl_vector_set( x , 0 , 0.0 );
		Stepper stepper;
		check_error_stepper_concept( stepper , constant_system2 , *x , *xerr );
		BOOST_CHECK_SMALL( fabs( gsl_vector_get( x , 0 ) - 0.1 ) , eps );
		gsl_vector_free( x ); gsl_vector_free( xerr );
	}
};

template< class Stepper >
struct perform_error_stepper_test< Stepper , vector_space_type >
{
	void operator()( void ) const
	{
		vector_space_type x , xerr;
		x.m_x = 0.0;
		Stepper stepper;
		check_error_stepper_concept( stepper , constant_system3 , x , xerr );
		BOOST_CHECK_SMALL( fabs( x.m_x - 0.1 ) , eps );
	}
};

template< class Stepper >
struct perform_error_stepper_test< Stepper , array_type >
{
	void operator()( void )
	{
		array_type x , xerr;
		x[0] = 0.0;
		Stepper stepper;
		check_error_stepper_concept( stepper , constant_system4 , x , xerr );
		BOOST_CHECK_SMALL( fabs( x[0] - 0.1 ) , eps );
	}
};



template< class State > class error_stepper_methods : public mpl::vector<
	explicit_error_rk54_ck< State , double , typename algebra_dispatcher< State >::type >
> { };

typedef mpl::insert_range<	mpl::vector0<> , mpl::end< mpl::vector0<> >::type ,	error_stepper_methods< vector_type > >::type first_error_stepper_type;
typedef mpl::insert_range<	first_error_stepper_type , mpl::end< first_error_stepper_type >::type , error_stepper_methods< gsl_vector_type > >::type second_error_stepper_type;
typedef mpl::insert_range<	second_error_stepper_type , mpl::end< second_error_stepper_type >::type , error_stepper_methods< vector_space_type > >::type third_error_stepper_type;
typedef mpl::insert_range<	third_error_stepper_type , mpl::end< third_error_stepper_type >::type , error_stepper_methods< array_type > >::type all_error_stepper_methods;


BOOST_AUTO_TEST_CASE_TEMPLATE( error_stepper_test , Stepper , all_error_stepper_methods )
{
	perform_error_stepper_test< Stepper , typename Stepper::state_type > tester;
	tester();
}
















/* ToDO: check actual results of controlled step... */


template< class ControlledStepper , class State > struct perform_controlled_stepper_test;

template< class ControlledStepper >
struct perform_controlled_stepper_test< ControlledStepper , vector_type >
{
	void operator()( void )
	{
		vector_type x( 1 , 0.0 );
		typename ControlledStepper::error_stepper_type error_stepper;
		error_checker_standard< typename ControlledStepper::state_type , typename ControlledStepper::time_type > error_checker;
		ControlledStepper controlled_stepper( error_stepper , error_checker );
		check_controlled_stepper_concept( controlled_stepper , constant_system1 , x );
		//BOOST_CHECK_SMALL( fabs( x[0] - 0.1 ) , eps );
	}
};

template< class ControlledStepper >
struct perform_controlled_stepper_test< ControlledStepper , gsl_vector_type >
{
	void operator()( void ) const
	{
		gsl_vector_type *x = gsl_vector_alloc( 1 );
		gsl_vector_set( x , 0 , 0.0 );
		typename ControlledStepper::error_stepper_type error_stepper;
		error_checker_standard< typename ControlledStepper::state_type , typename ControlledStepper::time_type > error_checker;
		ControlledStepper controlled_stepper( error_stepper , error_checker );
		check_controlled_stepper_concept( controlled_stepper , constant_system2 , *x );
		//BOOST_CHECK_SMALL( fabs( gsl_vector_get( x , 0 ) - 0.1 ) , eps );
		gsl_vector_free( x );
	}
};

template< class ControlledStepper >
struct perform_controlled_stepper_test< ControlledStepper , vector_space_type >
{
	void operator()( void ) const
	{
		vector_space_type x;
		x.m_x = 0.0;
		typename ControlledStepper::error_stepper_type error_stepper;
		error_checker_standard< typename ControlledStepper::state_type , typename ControlledStepper::time_type , vector_space_algebra > error_checker;
		ControlledStepper controlled_stepper( error_stepper , error_checker );
		check_controlled_stepper_concept( controlled_stepper , constant_system3 , x );
		//BOOST_CHECK_SMALL( fabs( x.m_x - 0.1 ) , eps );
	}
};

template< class ControlledStepper >
struct perform_controlled_stepper_test< ControlledStepper , array_type >
{
	void operator()( void )
	{
		array_type x;
		x[0] = 0.0;
		typename ControlledStepper::error_stepper_type error_stepper;
		error_checker_standard< typename ControlledStepper::state_type , typename ControlledStepper::time_type > error_checker;
		ControlledStepper controlled_stepper( error_stepper , error_checker );
		check_controlled_stepper_concept( controlled_stepper , constant_system4 , x );
		//BOOST_CHECK_SMALL( fabs( x[0] - 0.1 ) , eps );
	}
};


template< class State > class controlled_stepper_methods : public mpl::vector<
	controlled_error_stepper< explicit_error_rk54_ck< State , double , typename algebra_dispatcher< State >::type > >
> { };

typedef mpl::insert_range<	mpl::vector0<> , mpl::end< mpl::vector0<> >::type ,	controlled_stepper_methods< vector_type > >::type first_controlled_stepper_type;
typedef mpl::insert_range<	first_controlled_stepper_type , mpl::end< first_controlled_stepper_type >::type , controlled_stepper_methods< gsl_vector_type > >::type second_controlled_stepper_type;
//typedef mpl::insert_range<	second_controlled_stepper_type , mpl::end< second_controlled_stepper_type >::type , controlled_stepper_methods< array_type > >::type all_controlled_stepper_methods;
typedef mpl::insert_range<	second_controlled_stepper_type , mpl::end< second_controlled_stepper_type >::type , controlled_stepper_methods< vector_space_type > >::type third_controlled_stepper_type;
typedef mpl::insert_range<	third_controlled_stepper_type , mpl::end< third_controlled_stepper_type >::type , controlled_stepper_methods< array_type > >::type all_controlled_stepper_methods;

BOOST_AUTO_TEST_CASE_TEMPLATE( controlled_stepper_test , ControlledStepper , all_controlled_stepper_methods )
{
	perform_controlled_stepper_test< ControlledStepper , typename ControlledStepper::state_type > tester;
	tester();
}
