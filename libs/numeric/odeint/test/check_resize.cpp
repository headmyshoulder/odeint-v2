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
#include <boost/bind.hpp>
#include <boost/utility.hpp>

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

template< class Stepper >
void test_resize( Stepper &stepper , size_t multiplicity , size_t resize_calls )
{
	adjust_size_count = 0;

	test_array_type x;
	stepper.do_step( constant_system , x , 0.0 , 0.1 );
	stepper.do_step( constant_system , x , 0.0 , 0.1 );
	stepper.do_step( constant_system , x , 0.0 , 0.1 );

	BOOST_TEST_MESSAGE( "adjust_size_count : " << adjust_size_count );

	BOOST_CHECK_MESSAGE( adjust_size_count == resize_calls * multiplicity , "adjust_size_count : " << adjust_size_count );
}


BOOST_AUTO_TEST_SUITE( check_resize_test )

//test_suite* init_unit_test_suite( int argc, char* argv[] )
//{
//    test_suite *test = BOOST_TEST_SUITE( "check resize functionality" );

    typedef explicit_euler< test_array_type , double , standard_algebra< test_array_type > , standard_operations< double > , adjust_size_manually_tag > euler_manual_type;
    typedef explicit_euler< test_array_type , double , standard_algebra< test_array_type > , standard_operations< double > , adjust_size_initially_tag > euler_initially_type;
    typedef explicit_euler< test_array_type , double , standard_algebra< test_array_type > , standard_operations< double > , adjust_size_always_tag > euler_always_type;

    typedef explicit_rk4< test_array_type , double , standard_algebra< test_array_type > , standard_operations< double > , adjust_size_manually_tag > rk4_manual_type;
    typedef explicit_rk4< test_array_type , double , standard_algebra< test_array_type > , standard_operations< double > , adjust_size_initially_tag > rk4_initially_type;
    typedef explicit_rk4< test_array_type , double , standard_algebra< test_array_type > , standard_operations< double > , adjust_size_always_tag > rk4_always_type;


    euler_manual_type euler_manual;
    euler_initially_type euler_initially;
    euler_always_type euler_always;

    rk4_manual_type rk4_manual;
    rk4_initially_type rk4_initially;
    rk4_always_type rk4_always;

//    BOOST_AUTO_TEST_CASE( boost::bind( &test_resize< euler_manual_type > , boost::ref( euler_manual ) , 1 , 0 ) );



//    test->add( BOOST_TEST_CASE( boost::bind( &test_resize< euler_manual_type > , boost::ref( euler_manual ) , 1 , 0 ) ) );
//    test->add( BOOST_TEST_CASE( boost::bind( &test_resize< euler_initially_type > , boost::ref( euler_initially ) , 1 , 1 ) ) );
//    test->add( BOOST_TEST_CASE( boost::bind( &test_resize< euler_always_type > , boost::ref( euler_always ) , 1 , 3 ) ) );
//
//    test->add( BOOST_TEST_CASE( boost::bind( &test_resize< rk4_manual_type > , boost::ref( rk4_manual ) , 5 , 0 ) ) );
//    test->add( BOOST_TEST_CASE( boost::bind( &test_resize< rk4_initially_type > , boost::ref( rk4_initially ) , 5 , 1 ) ) );
//    test->add( BOOST_TEST_CASE( boost::bind( &test_resize< rk4_always_type > , boost::ref( rk4_always ) , 5 , 3 ) ) );

//    return test;
//}
BOOST_AUTO_TEST_SUITE_END()
