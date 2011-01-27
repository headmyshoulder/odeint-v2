/* Boost stepper_euler.cpp test file

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 This file tests the use of the euler stepper

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#define BOOST_TEST_MODULE odeint_resize

#include <vector>
#include <cmath>

#include <boost/array.hpp>
#include <boost/bind.hpp>
#include <boost/utility.hpp>

#include <boost/test/unit_test.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/at.hpp>

#include <boost/numeric/odeint/stepper/explicit_euler.hpp>
#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

namespace mpl = boost::mpl;

size_t adjust_size_count;

typedef boost::array< double , 1 > test_array_type;

namespace boost { namespace numeric { namespace odeint {


	template<>	struct is_resizeable< test_array_type >
	{
		struct type : public boost::true_type { };
		const static bool value = type::value;
	};


	template<> bool adjust_size( const test_array_type &x1 , test_array_type &x2 )
	{
		adjust_size_count++;
		return false;
	}

} } }



void constant_system( const test_array_type &x , test_array_type &dxdt , double t ) { dxdt[0] = 1.0; }


BOOST_AUTO_TEST_SUITE( check_resize_test )


typedef explicit_euler< test_array_type , double , test_array_type , double , range_algebra , default_operations , adjust_size_manually_tag > euler_manual_type;
typedef explicit_euler< test_array_type , double , test_array_type , double , range_algebra , default_operations , adjust_size_initially_tag > euler_initially_type;
typedef explicit_euler< test_array_type , double , test_array_type , double , range_algebra , default_operations , adjust_size_always_tag > euler_always_type;

typedef explicit_rk4< test_array_type , double , test_array_type , double , range_algebra , default_operations , adjust_size_manually_tag > rk4_manual_type;
typedef explicit_rk4< test_array_type , double , test_array_type , double , range_algebra , default_operations , adjust_size_initially_tag > rk4_initially_type;
typedef explicit_rk4< test_array_type , double , test_array_type , double , range_algebra , default_operations , adjust_size_always_tag > rk4_always_type;


typedef mpl::vector<
	mpl::vector< euler_manual_type , mpl::int_<1> , mpl::int_<0> > ,
	mpl::vector< euler_initially_type , mpl::int_<1> , mpl::int_<1> > ,
	mpl::vector< euler_always_type , mpl::int_<1> , mpl::int_<3> > ,
	mpl::vector< rk4_manual_type , mpl::int_<5> , mpl::int_<0> > ,
	mpl::vector< rk4_initially_type , mpl::int_<5> , mpl::int_<1> > ,
	mpl::vector< rk4_always_type , mpl::int_<5> , mpl::int_<3> >
>::type resize_check_types;


BOOST_AUTO_TEST_CASE_TEMPLATE( test_resize , T, resize_check_types )
{
	typedef typename mpl::at< T , mpl::int_< 0 > >::type stepper_type;
	const size_t resize_calls = mpl::at< T , mpl::int_< 1 > >::type::value;
	const size_t multiplicity = mpl::at< T , mpl::int_< 2 > >::type::value;
	adjust_size_count = 0;

	stepper_type stepper;
	test_array_type x;
	stepper.do_step( constant_system , x , 0.0 , 0.1 );
	stepper.do_step( constant_system , x , 0.0 , 0.1 );
	stepper.do_step( constant_system , x , 0.0 , 0.1 );

	BOOST_TEST_MESSAGE( "adjust_size_count : " << adjust_size_count );
	BOOST_CHECK_MESSAGE( adjust_size_count == resize_calls * multiplicity , "adjust_size_count : " << adjust_size_count );
}


BOOST_AUTO_TEST_SUITE_END()
