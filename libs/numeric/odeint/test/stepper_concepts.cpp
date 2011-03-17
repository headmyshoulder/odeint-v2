/* Boost stepper_euler.cpp test file

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 This file tests the use of the all different steppers with several state types:
 std::vector< double >
 vector_space_1d< double >  (see vector_space_1d.hpp)
 std::tr1::array< double , 1 >

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#define BOOST_TEST_MODULE odeint_stepper_concepts

#include <vector>
#include <cmath>
#include <iostream>

#include <tr1/array>

#include <boost/test/unit_test.hpp>

#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <boost/utility.hpp>
#include <boost/type_traits/add_reference.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/mpl/end.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/inserter.hpp>

#include <boost/numeric/odeint/stepper/explicit_euler.hpp>
#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>
#include <boost/numeric/odeint/stepper/explicit_error_rk54_ck.hpp>
#include <boost/numeric/odeint/stepper/explicit_error_dopri5.hpp>
#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

#include "vector_space_1d.hpp"


using std::vector;

using namespace boost::unit_test;
using namespace boost::numeric::odeint;
namespace mpl = boost::mpl;


typedef std::vector< double > vector_type;
typedef vector_space_1d< double > vector_space_type;
typedef std::tr1::array< double , 1 > array_type;


const double result = 2.2;

typedef mpl::vector< vector_type , vector_space_type , array_type >::type container_types;


template< class State > struct algebra_dispatcher { typedef range_algebra type; };
template<> struct algebra_dispatcher< vector_space_type > { typedef vector_space_algebra type; };
template<> struct algebra_dispatcher< double > { typedef vector_space_algebra type; };

struct constant_system_vector_class
{
	void operator()( const vector_type &x , vector_type &dxdt , double t ) const
	{
		dxdt[0] = 1.0;
	}
};

struct constant_system_vector_space_class
{
	void operator()( const vector_space_type &x , vector_space_type &dxdt , double t ) const
	{
		dxdt.m_x = 1.0;
	}
};

struct constant_system_array_class
{
	void operator()( const array_type &x , array_type &dxdt , double t ) const
	{
		dxdt[0] = 1.0;
	}
};

void constant_system_vector( const vector_type &x , vector_type &dxdt , double t ) { dxdt[0] = 1.0; }
void constant_system_vector_space( const vector_space_type &x , vector_space_type &dxdt , double t ) { dxdt.m_x = 1.0; }
void constant_system_array( const array_type &x , array_type &dxdt , double t ) { dxdt[0] = 1.0; }

const double eps = 1.0e-14;

template< class Stepper , class System >
void check_stepper_concept( Stepper &stepper , System system , typename Stepper::deriv_type &x )
{
    typedef Stepper stepper_type;
    typedef typename stepper_type::deriv_type container_type;
    typedef typename stepper_type::order_type order_type;
    typedef typename stepper_type::time_type time_type;

    stepper.do_step( system , x , 0.0 , 0.1 );
}

template< class Stepper , class System >
void check_error_stepper_concept( Stepper &stepper , System system , typename Stepper::state_type &x , typename Stepper::state_type &xerr )
{
    typedef Stepper stepper_type;
    typedef typename stepper_type::deriv_type container_type;
    typedef typename stepper_type::order_type order_type;
    typedef typename stepper_type::time_type time_type;

    stepper.do_step( system , typename boost::add_reference< container_type>::type( x ), 0.0 , 0.1 ,  typename boost::add_reference< container_type>::type( xerr ) );
}

template< class Stepper , class System >
void check_controlled_stepper_concept( Stepper &stepper , System system , typename Stepper::state_type &x )
{
	typedef Stepper stepper_type;
    typedef typename stepper_type::deriv_type container_type;
    typedef typename stepper_type::order_type order_type;
    typedef typename stepper_type::time_type time_type;

    time_type t = 0.0 , dt = 0.1;
    controlled_step_result step_result = stepper.try_step( system , x , t , dt );
    BOOST_CHECK_MESSAGE( step_result == success_step_size_increased , "step result: " << step_result ); // error = 0 for constant system -> step size is always too small
}





template< class Stepper , class State > struct perform_stepper_test;

template< class Stepper >
struct perform_stepper_test< Stepper , vector_type >
{
	void operator()( void )
	{
		vector_type x( 1 , 2.0 );
		Stepper stepper;
		check_stepper_concept( stepper , constant_system_vector , x );
		check_stepper_concept( stepper , boost::cref( constant_system_vector_class() ) , x );
		BOOST_CHECK_SMALL( fabs( x[0] - result ) , eps );
	}
};

template< class Stepper >
struct perform_stepper_test< Stepper , vector_space_type >
{
	void operator()( void ) const
	{
		vector_space_type x;
		x.m_x = 2.0;
		Stepper stepper;
		check_stepper_concept( stepper , constant_system_vector_space , x );
		check_stepper_concept( stepper , boost::cref( constant_system_vector_space_class() ) , x );
		BOOST_CHECK_SMALL( fabs( x.m_x - result ) , eps );
	}
};

template< class Stepper >
struct perform_stepper_test< Stepper , array_type >
{
	void operator()( void )
	{
		array_type x;
		x[0] = 2.0;
		Stepper stepper;
		check_stepper_concept( stepper , constant_system_array , x );
		check_stepper_concept( stepper , boost::cref( constant_system_array_class() ) , x );
		BOOST_CHECK_SMALL( fabs( x[0] - result ) , eps );
	}
};


template< class State > class stepper_methods : public mpl::vector<
	explicit_euler< State , double , State , double , typename algebra_dispatcher< State >::type > ,
	explicit_rk4< State , double , State , double , typename algebra_dispatcher< State >::type > ,
	explicit_error_rk54_ck< State , double , State , double , typename algebra_dispatcher< State >::type > ,
	explicit_error_dopri5< State , double , State , double , typename algebra_dispatcher< State >::type >
> { };



typedef mpl::copy
<
  container_types ,
  mpl::inserter
  <
    mpl::vector0<> ,
    mpl::insert_range
    <
      mpl::_1 ,
      mpl::end< mpl::_1 > ,
      stepper_methods< mpl::_2 >
    >
  >
>::type all_stepper_methods;



BOOST_AUTO_TEST_SUITE( stepper_concept_test )

BOOST_AUTO_TEST_CASE_TEMPLATE( stepper_test , Stepper, all_stepper_methods )
{
	perform_stepper_test< Stepper , typename Stepper::deriv_type > tester;
	tester();
}


BOOST_AUTO_TEST_SUITE_END()

















template< class Stepper , class State > struct perform_error_stepper_test;

template< class Stepper >
struct perform_error_stepper_test< Stepper , vector_type >
{
	void operator()( void )
	{
		vector_type x( 1 , 2.0 ) , xerr( 1 );
		Stepper stepper;
		check_error_stepper_concept( stepper , constant_system_vector , x , xerr );
		check_error_stepper_concept( stepper , boost::cref( constant_system_vector_class() ) , x , xerr );
		BOOST_CHECK_SMALL( fabs( x[0] - result ) , eps );
	}
};


template< class Stepper >
struct perform_error_stepper_test< Stepper , vector_space_type >
{
	void operator()( void ) const
	{
		vector_space_type x , xerr;
		x.m_x = 2.0;
		Stepper stepper;
		check_error_stepper_concept( stepper , constant_system_vector_space , x , xerr );
		check_error_stepper_concept( stepper , boost::cref( constant_system_vector_space_class() ) , x , xerr );
		BOOST_CHECK_SMALL( fabs( x.m_x - result ) , eps );
	}
};

template< class Stepper >
struct perform_error_stepper_test< Stepper , array_type >
{
	void operator()( void )
	{
		array_type x , xerr;
		x[0] = 2.0;
		Stepper stepper;
		check_error_stepper_concept( stepper , constant_system_array , x , xerr );
		check_error_stepper_concept( stepper , boost::cref( constant_system_array_class() ) , x , xerr );
		BOOST_CHECK_SMALL( fabs( x[0] - result ) , eps );
	}
};


template< class State > class error_stepper_methods : public mpl::vector<
	explicit_error_rk54_ck< State , double , State , double , typename algebra_dispatcher< State >::type > ,
	explicit_error_dopri5< State , double , State , double , typename algebra_dispatcher< State >::type >
> { };


typedef mpl::copy
<
  container_types ,
  mpl::inserter
  <
    mpl::vector0<> ,
    mpl::insert_range
    <
      mpl::_1 ,
      mpl::end< mpl::_1 > ,
      error_stepper_methods< mpl::_2 >
    >
  >
>::type all_error_stepper_methods;


BOOST_AUTO_TEST_SUITE( error_stepper_concept_test )

BOOST_AUTO_TEST_CASE_TEMPLATE( error_stepper_test , Stepper , all_error_stepper_methods )
{
	perform_error_stepper_test< Stepper , typename Stepper::state_type > tester;
	tester();
}

BOOST_AUTO_TEST_SUITE_END()
















/* ToDO: check actual results of controlled step... */


template< class ControlledStepper , class State > struct perform_controlled_stepper_test;

template< class ControlledStepper >
struct perform_controlled_stepper_test< ControlledStepper , vector_type >
{
	void operator()( void )
	{
		vector_type x( 1 , 2.0 );
		typename ControlledStepper::stepper_type error_stepper;
		default_error_checker< typename ControlledStepper::value_type > error_checker;
		ControlledStepper controlled_stepper( error_stepper , error_checker );
		check_controlled_stepper_concept( controlled_stepper , constant_system_vector , x );
		check_controlled_stepper_concept( controlled_stepper , boost::cref( constant_system_vector_class() ) , x );
		BOOST_CHECK_SMALL( fabs( x[0] - result ) , eps );
	}
};

template< class ControlledStepper >
struct perform_controlled_stepper_test< ControlledStepper , vector_space_type >
{
	void operator()( void ) const
	{
		vector_space_type x;
		x.m_x = 2.0;
		typename ControlledStepper::stepper_type error_stepper;
		default_error_checker< typename ControlledStepper::value_type , vector_space_algebra > error_checker;
		ControlledStepper controlled_stepper( error_stepper , error_checker );
		check_controlled_stepper_concept( controlled_stepper , constant_system_vector_space , x );
		check_controlled_stepper_concept( controlled_stepper , boost::cref( constant_system_vector_space_class() ) , x );
		BOOST_CHECK_SMALL( fabs( x.m_x - result ) , eps );
	}
};

template< class ControlledStepper >
struct perform_controlled_stepper_test< ControlledStepper , array_type >
{
	void operator()( void )
	{
		array_type x;
		x[0] = 2.0;
		typename ControlledStepper::stepper_type error_stepper;
		default_error_checker< typename ControlledStepper::value_type > error_checker;
		ControlledStepper controlled_stepper( error_stepper , error_checker );
		check_controlled_stepper_concept( controlled_stepper , constant_system_array , x );
		check_controlled_stepper_concept( controlled_stepper , boost::cref( constant_system_array_class() ) , x );
		BOOST_CHECK_SMALL( fabs( x[0] - result ) , eps );
	}
};

template< class State > class controlled_stepper_methods : public mpl::vector<
	controlled_error_stepper< explicit_error_rk54_ck< State , double , State , double , typename algebra_dispatcher< State >::type > > ,
	controlled_error_stepper< explicit_error_dopri5< State , double , State , double , typename algebra_dispatcher< State >::type > >
> { };

typedef mpl::copy
<
  container_types ,
  mpl::inserter
  <
    mpl::vector0<> ,
    mpl::insert_range
    <
      mpl::_1 ,
      mpl::end< mpl::_1 > ,
      controlled_stepper_methods< mpl::_2 >
    >
  >
>::type all_controlled_stepper_methods;




BOOST_AUTO_TEST_SUITE( controlled_stepper_concept_test )

BOOST_AUTO_TEST_CASE_TEMPLATE( controlled_stepper_test , ControlledStepper , all_controlled_stepper_methods )
{
	perform_controlled_stepper_test< ControlledStepper , typename ControlledStepper::state_type > tester;
	tester();
}

BOOST_AUTO_TEST_SUITE_END()
