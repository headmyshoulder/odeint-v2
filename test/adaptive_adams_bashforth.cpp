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

BOOST_AUTO_TEST_CASE( test_instantiation )
{
	typedef double time_type;
	typedef std::vector<double> state_type;

	adaptive_adams_bashforth<1, state_type> s1;
	adaptive_adams_bashforth<2, state_type> s2;
	adaptive_adams_bashforth<3, state_type> s3;
	adaptive_adams_bashforth<4, state_type> s4;
	adaptive_adams_bashforth<5, state_type> s5;
	adaptive_adams_bashforth<6, state_type> s6;
	adaptive_adams_bashforth<7, state_type> s7;
	adaptive_adams_bashforth<8, state_type> s8;
	adaptive_adams_bashforth<9, state_type> s9;

	state_type x0;
	x0.push_back(0);
	time_type t0 = 0.0;
	time_type dt = 0.1;

	s1.do_step(const_sys(), x0, t0, dt);
	s2.do_step(const_sys(), x0, t0, dt);
	s3.do_step(const_sys(), x0, t0, dt);
	s4.do_step(const_sys(), x0, t0, dt);
	s5.do_step(const_sys(), x0, t0, dt);
	s6.do_step(const_sys(), x0, t0, dt);
	s7.do_step(const_sys(), x0, t0, dt);
	s8.do_step(const_sys(), x0, t0, dt);
	s9.do_step(const_sys(), x0, t0, dt);
}

BOOST_AUTO_TEST_SUITE_END()