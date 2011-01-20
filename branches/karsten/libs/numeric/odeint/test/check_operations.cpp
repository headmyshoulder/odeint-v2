/*
 * check_operations.cpp
 *
 *  Created on: Jan 20, 2011
 *      Author: karsten
 */


#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <boost/numeric/odeint/algebra/standard_operations.hpp>

using boost::numeric::odeint::standard_operations;

const double eps = 1.0e-13;

struct double_fixture
{
	double_fixture( void )
	: res( 0.0 ) , x1( 1.0 ) , x2( 2.0 ) , x3( 3.0 ) , x4( 4.0 ) , x5( 5.0 ) , x6( 6.0 ) , x7( 7.0 ) , x8( 8.0 )
	{}
	~double_fixture( void )
	{
		BOOST_CHECK_CLOSE( x1 , double( 1.0 ) , eps );
		BOOST_CHECK_CLOSE( x2 , double( 2.0 ) , eps );
		BOOST_CHECK_CLOSE( x3 , double( 3.0 ) , eps );
		BOOST_CHECK_CLOSE( x4 , double( 4.0 ) , eps );
		BOOST_CHECK_CLOSE( x5 , double( 5.0 ) , eps );
		BOOST_CHECK_CLOSE( x6 , double( 6.0 ) , eps );
		BOOST_CHECK_CLOSE( x7 , double( 7.0 ) , eps );
		BOOST_CHECK_CLOSE( x8 , double( 8.0 ) , eps );
	}
	double res;
	double x1 , x2 , x3 , x4 , x5 , x6 , x7 , x8;
};

struct unit_fixture
{
	unit_fixture( void )
	{}
	~unit_fixture( void )
	{

	}


};


BOOST_AUTO_TEST_SUITE( check_operations_test )

BOOST_AUTO_TEST_CASE( scale_sum2_with_double )
{
	double_fixture f;
	double res = 0.0;
	typedef standard_operations::scale_sum2< double , double > Op;
	Op op( 1.25 , 9.81 );
	op( res , f.x1 , f.x2 );
	BOOST_CHECK_CLOSE( res , double( 20.87 ) , eps );
}

BOOST_AUTO_TEST_CASE( scale_sum3_with_double )
{
	double_fixture f;
	double res = 0.0;
	typedef standard_operations::scale_sum3< double , double , double > Op;
	Op op( 1.25 , 9.81 , 0.87 );
	op( res , f.x1 , f.x2 , f.x3 );
	BOOST_CHECK_CLOSE( res , double( 23.48 ) , eps );
}

BOOST_AUTO_TEST_CASE( scale_sum4_with_double )
{
	double_fixture f;
	double res = 0.0;
	typedef standard_operations::scale_sum4< double , double , double , double > Op;
	Op op( 1.25 , 9.81 , 0.87 , -0.15 );
	op( res , f.x1 , f.x2 , f.x3 , f.x4 );
	BOOST_CHECK_CLOSE( res , double( 22.88 ) , eps );
}

BOOST_AUTO_TEST_CASE( scale_sum5_with_double )
{
	double_fixture f;
	double res = 0.0;
	typedef standard_operations::scale_sum5< double , double , double , double , double > Op;
	Op op( 1.25 , 9.81 , 0.87 , -0.15 , -3.3 );
	op( res , f.x1 , f.x2 , f.x3 , f.x4 , f.x5 );
	BOOST_CHECK_CLOSE( res , double( 6.38 ) , eps );
}

BOOST_AUTO_TEST_CASE( scale_sum6_with_double )
{
	double_fixture f;
	double res = 0.0;
	typedef standard_operations::scale_sum6< double , double , double , double , double , double > Op;
	Op op( 1.25 , 9.81 , 0.87 , -0.15 , -3.3 , 4.2 );
	op( res , f.x1 , f.x2 , f.x3 , f.x4 , f.x5 , f.x6 );
	BOOST_CHECK_CLOSE( res , double( 31.58 ) , eps );
}

BOOST_AUTO_TEST_CASE( scale_sum7_with_double )
{
	double_fixture f;
	double res = 0.0;
	typedef standard_operations::scale_sum7< double , double , double , double , double , double , double > Op;
	Op op( 1.25 , 9.81 , 0.87 , -0.15 , -3.3 , 4.2 , -0.22 );
	op( res , f.x1 , f.x2 , f.x3 , f.x4 , f.x5 , f.x6 , f.x7 );
	BOOST_CHECK_CLOSE( res , double( 30.04 ) , eps );
}

BOOST_AUTO_TEST_CASE( rel_error_with_double )
{
	double_fixture f;
	double res = -1.1;
	typedef standard_operations::rel_error< double > Op;
	Op op( 0.1 , 0.2 , 0.15 , 0.12 );
	op( -f.x1 , -f.x2 , res );
	BOOST_CHECK_CLOSE( res , 6.17978 , 1.0e-4 );
}

BOOST_AUTO_TEST_CASE( maximum )
{
	double_fixture f;
	double res = 0.0;
	typedef standard_operations::maximum< double > Op;
	Op op;
	res = op( f.x1 , f.x2 );
	BOOST_CHECK_CLOSE( res , 2.0 , eps );
}

BOOST_AUTO_TEST_CASE( scale_sum2_with_units )
{

}



BOOST_AUTO_TEST_SUITE_END()
