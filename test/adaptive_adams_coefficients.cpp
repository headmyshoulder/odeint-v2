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

    typedef detail::adaptive_adams_coefficients<steps, deriv_type, time_type> aac_type;

	std::vector<double> deriv;
	deriv.push_back(-1);

	aac_type coeff;
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

BOOST_AUTO_TEST_CASE( test_copy )
{
	typedef std::vector<double> deriv_type;
	typedef double time_type;

	typedef detail::adaptive_adams_coefficients<3, deriv_type, time_type> aac_type;
	aac_type c1;

	deriv_type deriv(1);
	deriv[0] = 1.0;

	c1.step(deriv, 0.0);
	c1.confirm();
	c1.step(deriv, 1.0);
	c1.confirm();
	c1.step(deriv, 2.0);
	c1.confirm();

	aac_type c2(c1);
	BOOST_CHECK_EQUAL(c1.m_ss[0][0].m_v[0], c2.m_ss[0][0].m_v[0]);
	BOOST_CHECK(&(c1.m_ss[0][0].m_v) != &(c2.m_ss[0][0].m_v));

	aac_type c3;
	deriv_type *p1 = &(c3.m_ss[0][0].m_v);

	c3 = c1;
	// BOOST_CHECK(p1 == (&(c3.m_ss[0][0].m_v)));
	BOOST_CHECK_EQUAL(c1.m_ss[0][0].m_v[0], c3.m_ss[0][0].m_v[0]);
}

BOOST_AUTO_TEST_SUITE_END()