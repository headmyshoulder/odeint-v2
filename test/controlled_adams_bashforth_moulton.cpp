#include <boost/config.hpp>
#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif

#define BOOST_TEST_MODULE odeint_controlled_adams_bashforth_moulton

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/stepper/controlled_adams_bashforth_moulton.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

BOOST_AUTO_TEST_SUITE( controlled_adams_bashforth_moulton_test )

BOOST_AUTO_TEST_CASE( test_init )
{
	
}

BOOST_AUTO_TEST_SUITE_END()