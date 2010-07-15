/* Boost stepper_euler.cpp test file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file tests the use of the euler stepper
  
 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <vector>
#include <cmath>
#include <boost/array.hpp>

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

size_t adjust_size_count;

typedef boost::array< double , 1 > test_array_type;

namespace boost {
namespace numeric {
namespace odeint {



template<>
struct is_resizeable< test_array_type >
{
	struct type : public boost::true_type { };
	const static bool value = type::value;
};


template<>
void adjust_size( const test_array_type &x1 , test_array_type &x2 )
{
	adjust_size_count++;
}

}
}
}



void constant_system( const test_array_type &x , test_array_type &dxdt , double t )
{
	dxdt[0] = 1.0;
}


void test_manual_resize( void )
{
	adjust_size_count = 0;

	test_array_type x;
	explicit_euler< test_array_type , double , standard_algebra< test_array_type > , standard_operations< double > , adjust_size_manually_tag > euler;
	euler.do_step( constant_system , x , 0.0 , 0.1 );

	BOOST_CHECK( adjust_size_count == 0 );
}

void test_initially_resize( void )
{
	adjust_size_count = 0;
	test_array_type x;
	explicit_euler< test_array_type , double , standard_algebra< test_array_type > , standard_operations< double > , adjust_size_initially_tag > euler;
	euler.do_step( constant_system , x , 0.0 , 0.1 );
	euler.do_step( constant_system , x , 0.0 , 0.1 );
	euler.do_step( constant_system , x , 0.0 , 0.1 );
	BOOST_CHECK( adjust_size_count == 1 );
}

void test_always_resize( void )
{
	adjust_size_count = 0;
	test_array_type x;
	explicit_euler< test_array_type , double , standard_algebra< test_array_type > , standard_operations< double > , adjust_size_always_tag > euler;
	euler.do_step( constant_system , x , 0.0 , 0.1 );
	euler.do_step( constant_system , x , 0.0 , 0.1 );
	euler.do_step( constant_system , x , 0.0 , 0.1 );
	BOOST_CHECK( adjust_size_count == 3 );
}



test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE( "check resize functionality" );

    test->add( BOOST_TEST_CASE( &test_manual_resize ) );
    test->add( BOOST_TEST_CASE( &test_initially_resize ) );
    test->add( BOOST_TEST_CASE( &test_always_resize ) );

    return test;
}

