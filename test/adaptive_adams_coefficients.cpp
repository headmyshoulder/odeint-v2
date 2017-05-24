#include <boost/config.hpp>
#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif

#define BOOST_TEST_MODULE odeint_adaptive_adams_coefficients

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/stepper/detail/adaptive_adams_coefficients.hpp>

#include <vector>

#include <boost/mpl/list.hpp>
#include <boost/mpl/size_t.hpp>
#include <boost/mpl/range_c.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

typedef double value_type;

BOOST_AUTO_TEST_SUITE( adaptive_adams_coefficients_test )

typedef boost::mpl::range_c< size_t , 2 , 10 > vector_of_steps;
BOOST_AUTO_TEST_CASE_TEMPLATE( test_step, step_type, vector_of_steps )
{
	const static size_t steps = step_type::value;

	typedef std::vector<double> deriv_type;
	typedef double time_type;

    typedef detail::adaptive_adams_coefficients<steps, deriv_type, time_type> aac;

	std::vector<double> deriv;
	deriv.push_back(-1);

	aac coeff;
	for(size_t i=0; i<steps; ++i)
	{
		coeff.step(deriv, i);
		BOOST_CHECK_EQUAL(coeff.m_tts[0], i);

		coeff.confirm();
		BOOST_CHECK_EQUAL(coeff.m_ts[0], i);
	}

	BOOST_CHECK_EQUAL(coeff.m_ss[0][0].m_v[0], -1);
	BOOST_CHECK_EQUAL(coeff.m_ss[1][0].m_v[0], 0);
}

BOOST_AUTO_TEST_SUITE_END()