/*
 [auto_generated]
 libs/numeric/odeint/test/symplectic_steppers.cpp

 [begin_description]
 tba.
 [end_description]

 Copyright 2009-2012 Karsten Ahnert
 Copyright 2009-2012 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/config.hpp>
#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif

#define BOOST_TEST_MODULE odeint_symplectic_steppers

#include <boost/test/unit_test.hpp>

#include <boost/array.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/zip_view.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/mpl/end.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/inserter.hpp>
#include <boost/mpl/at.hpp>
// #include <boost/mpl/for_each.hpp>
#include <boost/fusion/container/vector.hpp>

#include <boost/numeric/odeint/stepper/symplectic_euler.hpp>
#include <boost/numeric/odeint/stepper/symplectic_rkn_sb3a_mclachlan.hpp>
#include <boost/numeric/odeint/stepper/symplectic_rkn_sb3a_m4_mclachlan.hpp>

#include "diagnostic_state_type.hpp"
#include "const_range.hpp"

#include <iostream>
#include <typeinfo>
using namespace std;



using namespace boost::unit_test;
using namespace boost::numeric::odeint;
namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

class custom_range_algebra : public range_algebra { };
class custom_default_operations : public default_operations { };


template< class Resizer >
class vector_steppers : public mpl::vector<
    symplectic_euler< diagnostic_state_type , diagnostic_state_type2 , double ,
                      diagnostic_deriv_type , diagnostic_deriv_type2 , double ,
                      custom_range_algebra , custom_default_operations , Resizer > ,
    symplectic_rkn_sb3a_mclachlan< diagnostic_state_type , diagnostic_state_type2 , double ,
                                   diagnostic_deriv_type , diagnostic_deriv_type2 , double ,
                                   custom_range_algebra , custom_default_operations , Resizer > ,
    symplectic_rkn_sb3a_m4_mclachlan< diagnostic_state_type , diagnostic_state_type2 , double ,
                                      diagnostic_deriv_type , diagnostic_deriv_type2 , double ,
                                      custom_range_algebra , custom_default_operations , Resizer > 
    > {};

typedef mpl::vector< initially_resizer , always_resizer , never_resizer > resizers;
typedef mpl::vector_c< int , 1 , 3 , 0 > resizer_multiplicities ;


typedef mpl::copy<
    resizers ,
    mpl::inserter<
        mpl::vector0<> ,
        mpl::insert_range<
            mpl::_1 ,
            mpl::end< mpl::_1 > ,
            vector_steppers< mpl::_2 >
            >
        >
    >::type all_stepper_methods;


typedef mpl::size< vector_steppers< initially_resizer > >::type num_steppers;
typedef mpl::copy<
    resizer_multiplicities ,
    mpl::inserter<
        mpl::vector0<> ,
        mpl::insert_range<
            mpl::_1 , 
            mpl::end< mpl::_1 > ,
            const_range< num_steppers , mpl::_2 >
            >
        >
    >::type all_multiplicities;
                                                                                 

struct constant_system
{
    template< class StateIn , class StateOut >
    void operator()( const StateIn &q , StateOut &dp ) const
    {
        dp[0] = 1.0;
    }
};




BOOST_AUTO_TEST_SUITE( symplectic_steppers_test )


BOOST_AUTO_TEST_CASE_TEMPLATE( test_assoc_types , Stepper , vector_steppers< initially_resizer > )
{
    BOOST_STATIC_ASSERT_MSG(
        ( boost::is_same< typename Stepper::coor_type , diagnostic_state_type >::value ) ,
        "Coordinate type" );
    BOOST_STATIC_ASSERT_MSG(
        ( boost::is_same< typename Stepper::momentum_type , diagnostic_state_type2 >::value ) ,
        "Momentum type" );
    BOOST_STATIC_ASSERT_MSG(
        ( boost::is_same< typename Stepper::coor_deriv_type , diagnostic_deriv_type >::value ) ,
        "Coordinate deriv type" );
    BOOST_STATIC_ASSERT_MSG(
        ( boost::is_same< typename Stepper::momentum_deriv_type , diagnostic_deriv_type2 >::value ) ,
        "Momentum deriv type" );

    BOOST_STATIC_ASSERT_MSG(
        ( boost::is_same< typename Stepper::state_type , std::pair< diagnostic_state_type , diagnostic_state_type2 > >::value ) ,
        "State type" );
    BOOST_STATIC_ASSERT_MSG(
        ( boost::is_same< typename Stepper::deriv_type , std::pair< diagnostic_deriv_type , diagnostic_deriv_type2 > >::value ) ,
        "Deriv type" );

    BOOST_STATIC_ASSERT_MSG( ( boost::is_same< typename Stepper::value_type , double >::value ) , "Value type" );
    BOOST_STATIC_ASSERT_MSG( ( boost::is_same< typename Stepper::time_type , double >::value ) , "Time type" );
    BOOST_STATIC_ASSERT_MSG( ( boost::is_same< typename Stepper::algebra_type , custom_range_algebra >::value ) , "Algebra type" );
    BOOST_STATIC_ASSERT_MSG( ( boost::is_same< typename Stepper::operations_type , custom_default_operations >::value ) , "Operations type" );

    BOOST_STATIC_ASSERT_MSG( ( boost::is_same< typename Stepper::resizer_type , initially_resizer >::value ) , "Resizer type" );
    BOOST_STATIC_ASSERT_MSG( ( boost::is_same< typename Stepper::stepper_category , stepper_tag >::value ) , "Stepper category" );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( test_adjust_size , Stepper , vector_steppers< initially_resizer > )
{
    counter_state::init_counter();
    counter_deriv::init_counter();
    counter_state2::init_counter();
    counter_deriv2::init_counter();

    {
        Stepper stepper;
        diagnostic_state_type x;
        stepper.adjust_size( x );
    }

    TEST_COUNTERS( counter_state , 0 , 0 , 0 , 0 );
    TEST_COUNTERS( counter_state2 , 0 , 0 , 0 , 0 );
    TEST_COUNTERS( counter_deriv , 1 , 1 , 0 , 1 );
    TEST_COUNTERS( counter_deriv2 , 1 , 1 , 0 , 1 );
}


typedef mpl::zip_view< mpl::vector< all_stepper_methods , all_multiplicities > >::type zipped_steppers;
BOOST_AUTO_TEST_CASE_TEMPLATE( test_resizing , Stepper , zipped_steppers )
{
    typedef typename mpl::at_c< Stepper , 0 >::type stepper_type;
    const size_t multiplicity = mpl::at_c< Stepper , 1 >::type::value;

    counter_state::init_counter();
    counter_deriv::init_counter();
    counter_state2::init_counter();
    counter_deriv2::init_counter();

    {
        stepper_type stepper;
        std::pair< diagnostic_state_type , diagnostic_state_type2 > x;
        stepper.do_step( constant_system() , x , 0.0 , 0.1 );
        stepper.do_step( constant_system() , x , 0.0 , 0.1 );
        stepper.do_step( constant_system() , x , 0.0 , 0.1 );
    }

    TEST_COUNTERS( counter_state , 0 , 0 , 0 , 0 );
    TEST_COUNTERS( counter_state2 , 0 , 0 , 0 , 0 );
    TEST_COUNTERS( counter_deriv , multiplicity , 1 , 0 , 1 );
    TEST_COUNTERS( counter_deriv2 , multiplicity , 1 , 0 , 1 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( test_copying1 , Stepper , vector_steppers< initially_resizer > )
{
    counter_state::init_counter();
    counter_deriv::init_counter();
    counter_state2::init_counter();
    counter_deriv2::init_counter();

    {
        Stepper stepper;
        Stepper stepper2( stepper );
    }

    TEST_COUNTERS( counter_state , 0 , 0 , 0 , 0 );
    TEST_COUNTERS( counter_state2 , 0 , 0 , 0 , 0 );
    TEST_COUNTERS( counter_deriv , 0 , 2 , 1 , 2 );
    TEST_COUNTERS( counter_deriv2 , 0 , 2 , 1 , 2 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( test_copying2 , Stepper , vector_steppers< initially_resizer > )
{
    counter_state::init_counter();
    counter_deriv::init_counter();
    counter_state2::init_counter();
    counter_deriv2::init_counter();

    {
        Stepper stepper;
        std::pair< diagnostic_state_type , diagnostic_state_type2 > x;
        stepper.do_step( constant_system() , x , 0.0 , 0.1 );
        Stepper stepper2( stepper );
    }

    TEST_COUNTERS( counter_state , 0 , 0 , 0 , 0 );
    TEST_COUNTERS( counter_state2 , 0 , 0 , 0 , 0 );
    TEST_COUNTERS( counter_deriv , 1 , 2 , 1 , 2 );
    TEST_COUNTERS( counter_deriv2 , 1 , 2 , 1 , 2 );
}




//[ put in separate file
// also test boost::ref
BOOST_AUTO_TEST_CASE( test_do_step_v1 )
{
}

BOOST_AUTO_TEST_CASE( test_do_step_v2 )
{
}
//]


BOOST_AUTO_TEST_CASE( test_with_vector_space_algebra )
{
}

BOOST_AUTO_TEST_CASE( test_with_fusion_algebra )
{
}

BOOST_AUTO_TEST_CASE( test_with_boost_units )
{
}





BOOST_AUTO_TEST_SUITE_END()
