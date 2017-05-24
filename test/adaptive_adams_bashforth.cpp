#include <boost/config.hpp>
#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif

#define BOOST_TEST_MODULE odeint_adaptive_adams_bashforth

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/stepper/adaptive_adams_bashforth.hpp>

#include <vector>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

struct const_sys
{
    template< class State , class Deriv , class Value >
    void operator()( const State &x , Deriv &dxdt , const Value &dt ) const
    {
        dxdt[0] = 1;
    }
};

BOOST_AUTO_TEST_SUITE( adaptive_adams_bashforth_test )

BOOST_AUTO_TEST_CASE( test_const_int )
{
	typedef double time_type;
	typedef std::vector<double> state_type;

	typedef adaptive_adams_bashforth<3, state_type, time_type> aab_type;

	state_type x0;
	x0.push_back(0);
	time_type t0 = 0.0;
	time_type dt = 0.1;

	aab_type aab;
	aab.do_step(const_sys(), x0, t0, dt);
}

BOOST_AUTO_TEST_SUITE_END()