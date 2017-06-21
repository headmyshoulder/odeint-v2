#include <boost/config.hpp>
#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif

#define BOOST_TEST_MODULE odeint_polynomial

#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/stepper/detail/polynomial.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint::detail;

BOOST_AUTO_TEST_SUITE( polynomial_test )

BOOST_AUTO_TEST_CASE( test_init )
{
	Polynomial<3, double> p1;

	boost::array<double, 3> coeff = {0, 1, 2};
	Polynomial<3, double> p2(coeff);	
}
BOOST_AUTO_TEST_CASE( test_add_roots )
{
	Polynomial<3, double> poly;
	poly.add_root(1);
	poly.add_root(0);
}
BOOST_AUTO_TEST_CASE( test_copy )
{
	typedef Polynomial<5, double> poly_type;
	poly_type p1;
	p1.add_root(1);
	p1.add_root(0);

	poly_type p2(p1);
	BOOST_CHECK_EQUAL(p1.m_coeff[0], p2.m_coeff[0]);
	BOOST_CHECK_EQUAL(p1.m_coeff[1], p2.m_coeff[1]);

	poly_type p3;
	double* a1 = &(p3.m_coeff[0]);

	p3 = p1;

	BOOST_CHECK(a1 == &(p3.m_coeff[0]));

	BOOST_CHECK_EQUAL(p1.m_coeff[0], p3.m_coeff[0]);
	BOOST_CHECK_EQUAL(p1.m_coeff[1], p3.m_coeff[1]);
}
BOOST_AUTO_TEST_CASE( test_remove )
{
	
}
BOOST_AUTO_TEST_CASE( test_integrate )
{
	
}
BOOST_AUTO_TEST_CASE( test_evaluate )
{
	boost::array<double, 3> c1 = {2, 1, 0};
	Polynomial<3, double> p1(c1);

	BOOST_CHECK_EQUAL(p1.evaluate(0), 0);
	BOOST_CHECK_EQUAL(p1.evaluate(1), 3);
	BOOST_CHECK_EQUAL(p1.evaluate(2), 10);

	boost::array<double, 5> c2 = {0.001, 10, 0, 3, 1};
	Polynomial<5, double> p2(c2);

	BOOST_CHECK_CLOSE(p2.evaluate(0), 1, 0.001);
	BOOST_CHECK_CLOSE(p2.evaluate(0.001), 1.003, 0.001);
	BOOST_CHECK_CLOSE(p2.evaluate(1.001), 14.034, 0.001);	
}

BOOST_AUTO_TEST_SUITE_END()