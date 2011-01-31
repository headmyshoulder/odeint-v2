/*
 * stepper_copying.cpp
 *
 *  Created on: Jan 23, 2011
 *      Author: karsten
 */


#define BOOST_TEST_MODULE odeint_stepper_copying

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/util/construct.hpp>
#include <boost/numeric/odeint/util/destruct.hpp>
#include <boost/numeric/odeint/util/copy.hpp>

#include <boost/numeric/odeint/stepper/explicit_euler.hpp>
#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>
#include <boost/numeric/odeint/stepper/explicit_error_rk54_ck.hpp>
#include <boost/numeric/odeint/stepper/explicit_error_dopri5.hpp>
#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>
#include <boost/numeric/odeint/stepper/dense_output_explicit.hpp>
#include <boost/numeric/odeint/stepper/dense_output_controlled_explicit_fsal.hpp>

template< class T , size_t Dim >
class test_array
{
public:

	const static size_t dim = Dim;
	typedef T value_type;
	typedef value_type* iterator;
	typedef const value_type* const_iterator;

	value_type& operator[]( size_t i )
	{
		return m_data[i];
	}

	const value_type& operator[]( size_t i ) const
	{
		return m_data[i];
	}

	iterator begin( void )
	{
		return m_data;
	}

	iterator end( void )
	{
		return m_data + dim;
	}

	const_iterator begin( void ) const
	{
		return m_data;
	}

	const_iterator end( void ) const
	{
		return m_data + dim;
	}


private:

	value_type m_data[dim];
};

template< class T , size_t Dim >
class test_array2 : public test_array< T , Dim >
{
};



/*
 * Explicit testing if copying was successful is difficult,
 * hence we only test if the number of copy operations is right.
 *
 * Otherwise one has to prepare the states.
 */

size_t construct_count = 0;
size_t construct2_count = 0;
size_t destruct_count = 0;
size_t destruct2_count = 0;
size_t copy_count = 0;
size_t copy2_count = 0;

void reset_counter( void )
{
	construct_count = 0;
	construct2_count = 0;
	destruct_count = 0;
	destruct2_count = 0;
	copy_count = 0;
	copy2_count = 0;
}


namespace boost { namespace numeric { namespace odeint {

template< class T , size_t Dim >
struct construct_impl< test_array< T , Dim > >
{
	static void construct( test_array< T , Dim > &arr )
	{
		++construct_count;
	}
};

template< class T , size_t Dim >
struct construct_impl< test_array2< T , Dim > >
{
	static void construct( test_array2< T , Dim > &arr )
	{
		++construct2_count;
	}
};


template< class T , size_t Dim >
struct destruct_impl< test_array< T , Dim > >
{
	static void destruct( test_array< T , Dim > &arr )
	{
		++destruct_count;
	}
};

template< class T , size_t Dim >
struct destruct_impl< test_array2< T , Dim > >
{
	static void destruct( test_array2< T , Dim > &arr )
	{
		++destruct2_count;
	}
};


template< class T , size_t Dim >
struct copy_impl< test_array< T , Dim > , test_array< T , Dim > >
{
	static void copy( const test_array< T , Dim > &from , test_array< T , Dim > &to )
	{
		++copy_count;
	}
};

template< class T , size_t Dim >
struct copy_impl< test_array2< T , Dim > , test_array2< T , Dim > >
{
	static void copy( const test_array2< T , Dim > &from , test_array2< T , Dim > &to )
	{
		++copy2_count;
	}
};


} } }



typedef test_array< double , 3 > state_type;
typedef test_array2< double , 3 > deriv_type;
typedef boost::numeric::odeint::explicit_euler< state_type , double , deriv_type > euler_type;
typedef boost::numeric::odeint::explicit_rk4< state_type , double , deriv_type > rk4_type;
typedef boost::numeric::odeint::explicit_error_rk54_ck< state_type , double , deriv_type > rk54_type;
typedef boost::numeric::odeint::explicit_error_dopri5< state_type , double , deriv_type > dopri5_type;
typedef boost::numeric::odeint::controlled_error_stepper< rk54_type > controlled_rk54_type;
typedef boost::numeric::odeint::controlled_error_stepper< dopri5_type > controlled_dopri5_type;
typedef boost::numeric::odeint::dense_output_explicit< euler_type > dense_output_euler_type;
typedef boost::numeric::odeint::dense_output_controlled_explicit_fsal< controlled_dopri5_type > dense_output_dopri5_type;

#define CHECK_COUNTERS( c1 , c2 , c3 , c4 , c5 , c6 ) \
	BOOST_CHECK_EQUAL( construct_count , size_t( c1 ) ); \
	BOOST_CHECK_EQUAL( construct2_count , size_t( c2 ) ); \
	BOOST_CHECK_EQUAL( destruct_count , size_t( c3 ) ); \
	BOOST_CHECK_EQUAL( destruct2_count , size_t( c4) ); \
	BOOST_CHECK_EQUAL( copy_count , size_t( c5 ) ) ; \
	BOOST_CHECK_EQUAL( copy2_count, size_t( c6 ) )

BOOST_AUTO_TEST_SUITE( stepper_copying )

/*
 * Construct + Destruct
 * 1 deriv_type in explicit_stepper_base
 */
BOOST_AUTO_TEST_CASE( explicit_euler_construct )
{
	reset_counter();
	{
		euler_type euler;
	}
	CHECK_COUNTERS( 0 , 1 , 0 , 1 , 0 , 0 );
}


/*
 * Construct + Destruct
 * 2 * 1 deriv_type in explicit_stepper_base
 *
 * Copying
 * 1 deriv_type in explicit_stepper_base
 */
BOOST_AUTO_TEST_CASE( explicit_euler_copy_construct )
{
	reset_counter();
	{
		euler_type euler;
		euler_type euler2( euler );
	}
	CHECK_COUNTERS( 0 , 1 + 1 , 0 , 1 + 1 , 0 , 1 );
}

/*
 * Construct + Destruct
 * 2 * 1 deriv_type in explicit_stepper_base
 *
 * Copying
 * 1 deriv_type in explicit_stepper_base
 */
BOOST_AUTO_TEST_CASE( explicit_euler_assign )
{
	reset_counter();
	{
		euler_type euler;
		euler_type euler2;
		euler2 = euler;
	}
	CHECK_COUNTERS( 0 , 2 , 0 , 2 , 0 , 1 );
}

/*
 * Construct + Destruct
 * 1 deriv_type in explicit_stepper_base
 * 3 deriv_type in explicit_rk4
 * 1 state_type in explicit_rk4
 */
BOOST_AUTO_TEST_CASE( explicit_rk4_construct )
{
	reset_counter();
	{
		rk4_type rk4;
	}
	CHECK_COUNTERS( 1 , 4 , 1 , 4 , 0 , 0 );
}

/*
 * Construct + Destruct
 * 2 * 1 deriv_type in explicit_stepper_base
 * 2 * 3 deriv_type in explicit_rk4
 * 2 * 1 state_type in explicit_rk4
 *
 * Copying
 * 1 deriv_type in explicit_stepper_base
 * 3 deriv_type in explicit_stepper_base
 * 1 state_type in explicit_stepper_base
 */
BOOST_AUTO_TEST_CASE( explicit_rk4_copy_construct )
{
	reset_counter();
	{
		rk4_type rk4;
		rk4_type rk4_2( rk4 );
	}
	CHECK_COUNTERS( 2 , 8 , 2 , 8 , 1 , 4 );
}

/*
 * Construct + Destruct
 * 2 * 1 deriv_type in explicit_stepper_base
 * 2 * 3 deriv_type in explicit_rk4
 * 2 * 1 state_type in explicit_rk4
 *
 * Copying
 * 1 deriv_type in explicit_stepper_base
 * 3 deriv_type in explicit_stepper_base
 * 1 state_type in explicit_stepper_base
 */
BOOST_AUTO_TEST_CASE( explicit_rk4_assign )
{
	reset_counter();
	{
		rk4_type rk4;
		rk4_type rk4_2;
		rk4 = rk4_2;
	}
	CHECK_COUNTERS( 2 , 8 , 2 , 8 , 1 , 4 );
}

/*
 * Construct + Destruct
 * 2 explicit_rk54_ck:
 * 2 * 1 deriv_type in explicit_error_stepper_base
 * 2 * 5 deriv_type in explicit_error_rk54_ck
 * 2 * 1 state_type in explicit_error_rk4
 * 1 controlled_stepper:
 * 1 deriv_type
 * 2 state_type
 *
 * Copying
 * 1 copy process of explicit_rk54_ck:
 * 1 deriv_type from explicit_error_stepper_base
 * 5 deriv_type from explicit_error_rk54_ck
 * 1 state_type from explicit_error_rk54_ck
 */
BOOST_AUTO_TEST_CASE( controlled_rk54_construct )
{
	reset_counter();
	{
		controlled_rk54_type stepper;
	}
	CHECK_COUNTERS( 4 , 13 , 4 , 13 , 1 , 6 );
}


/*
 * Construct + Destruct
 * 3 explicit_rk54_ck:
 * 3 * 1 deriv_type in explicit_error_stepper_base
 * 3 * 5 deriv_type in explicit_error_rk54_ck
 * 3 * 1 state_type in explicit_error_rk4
 * 2 controlled_stepper:
 * 2 * 1 deriv_type
 * 2 * 2 state_type
 *
 * Copying
 * 1 copy process of explicit_rk54_ck:
 * 1 deriv_type from explicit_error_stepper_base
 * 5 deriv_type from explicit_error_rk54_ck
 * 1 state_type from explicit_error_rk54_ck
 *
 * 1 process of copying controlled_error_stepper
 * 1 deriv_type from explicit_error_stepper_base
 * 5 deriv_type from explicit_error_rk54_ck
 * 1 state_type from explicit_error_rk54_ck
 * 1 deriv_type from controlled_error_stepper
 * 2 state_type from controlled_error_stepper
 */
BOOST_AUTO_TEST_CASE( controlled_rk54_copy_construct )
{
	reset_counter();
	{
		controlled_rk54_type stepper;
		controlled_rk54_type stepper2( stepper );
	}
	CHECK_COUNTERS( 7 , 20 , 7 , 20 , 4 , 13 );
}

/*
 * Construct + Destruct
 * 4 explicit_rk54_ck:
 * 4 * 1 deriv_type in explicit_error_stepper_base
 * 4 * 5 deriv_type in explicit_error_rk54_ck
 * 4 * 1 state_type in explicit_error_rk4
 * 2 controlled_stepper:
 * 2 * 1 deriv_type
 * 2 * 2 state_type
 *
 * Copying
 * 2 copy process of explicit_rk54_ck:
 * 2 * 1 deriv_type from explicit_error_stepper_base
 * 2 * 5 deriv_type from explicit_error_rk54_ck
 * 2 * 1 state_type from explicit_error_rk54_ck
 *
 * 1 process of copying controlled_error_stepper
 * 1 deriv_type from explicit_error_stepper_base
 * 5 deriv_type from explicit_error_rk54_ck
 * 1 state_type from explicit_error_rk54_ck
 * 1 deriv_type from controlled_error_stepper
 * 2 state_type from controlled_error_stepper
 */
BOOST_AUTO_TEST_CASE( controlled_rk54_assign )
{
	reset_counter();
	{
		controlled_rk54_type stepper;
		controlled_rk54_type stepper2;
		stepper2 = stepper;
	}
	CHECK_COUNTERS( 8 , 26 , 8 , 26 , 5 , 19 );
}


BOOST_AUTO_TEST_CASE( controlled_dopri5_construct )
{
	reset_counter();
	{
		controlled_dopri5_type dopri5;
	}
	// CHECK_COUNTERS( 1 , 1 , 1 , 1 , 1 );
}

BOOST_AUTO_TEST_CASE( controlled_dopri5_copy_construct )
{
	reset_counter();
	{
		controlled_dopri5_type dopri5;
		controlled_dopri5_type dopri5_2( dopri5 );
	}
	// CHECK_COUNTERS( 1 , 1 , 1 , 1 , 1 );
}

BOOST_AUTO_TEST_CASE( controlled_dopri5_assign )
{
	reset_counter();
	{
		controlled_dopri5_type dopri5;
		controlled_dopri5_type dopri5_2;
		dopri5_2 = dopri5;
	}
	// CHECK_COUNTERS( 1 , 1 , 1 , 1 , 1 );
}

BOOST_AUTO_TEST_CASE( dense_output_euler_construct )
{
	reset_counter();
	{
		dense_output_euler_type euler;
	}
	// CHECK_COUNTERS( 1 , 1 , 1 , 1 , 1 );
}

BOOST_AUTO_TEST_CASE( dense_output_euler_copy_construct )
{
	reset_counter();
	{
		dense_output_euler_type euler;
		dense_output_euler_type euler2( euler );
	}
	// CHECK_COUNTERS( 1 , 1 , 1 , 1 , 1 );
}

BOOST_AUTO_TEST_CASE( dense_output_euler_assign )
{
	reset_counter();
	{
		dense_output_euler_type euler;
		dense_output_euler_type euler2;
		euler2 = euler;
	}
	// CHECK_COUNTERS( 1 , 1 , 1 , 1 , 1 );
}

BOOST_AUTO_TEST_CASE( dense_output_dopri5_construct )
{
	reset_counter();
	{
		dense_output_dopri5_type dopri5;
	}
	// CHECK_COUNTERS( 1 , 1 , 1 , 1 , 1 );
}

BOOST_AUTO_TEST_CASE( dense_output_dopri5_copy_construct )
{
	reset_counter();
	{
		dense_output_dopri5_type dopri5;
		dense_output_dopri5_type dopri5_2( dopri5 );
	}
	// CHECK_COUNTERS( 1 , 1 , 1 , 1 , 1 );
}

BOOST_AUTO_TEST_CASE( dense_output_dopri5_assign )
{
	reset_counter();
	{
		dense_output_dopri5_type dopri5;
		dense_output_dopri5_type dopri5_2;
		dopri5_2 = dopri5;
	}
	// CHECK_COUNTERS( 1 , 1 , 1 , 1 , 1 );
}









BOOST_AUTO_TEST_SUITE_END()

