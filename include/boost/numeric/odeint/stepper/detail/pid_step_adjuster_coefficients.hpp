#ifndef PID_STEP_ADJUSTER_COEFFICIENTS_HPP_INCLUDED
#define PID_STEP_ADJUSTER_COEFFICIENTS_HPP_INCLUDED

#include <boost/array.hpp>

namespace boost{
namespace numeric{
namespace odeint{
namespace detail{

enum adjuster_type{
	BASIC = 1,
	H0211 = 2,
	H211b = 3,
	H211PI = 4,
	H0312 = 5,
	H312b = 6,
	H312PID = 7,
	H0321 = 8,
	H321 = 9
};

template<int Type >
class pid_step_adjuster_coefficients;

template<>
class pid_step_adjuster_coefficients<adjuster_type::BASIC> : public boost::array<double, 5>
{
public:
	pid_step_adjuster_coefficients()
	: boost::array<double, 5>({1, 0, 0, 0, 0})
	{}
};

template<>
class pid_step_adjuster_coefficients<adjuster_type::H0211> : public boost::array<double, 5>
{
public:
	pid_step_adjuster_coefficients()
	: boost::array<double, 5>({1.0/2, 1.0/2, 0, 1.0/2, 0})
	{}
};

template<>
class pid_step_adjuster_coefficients<adjuster_type::H211b> : public boost::array<double, 5>
{
public:
	pid_step_adjuster_coefficients()
	: boost::array<double, 5>({1.0/5, 2.0/5, 0, 1.0/5, 0})
	{}
};

template<>
class pid_step_adjuster_coefficients<adjuster_type::H211PI> : public boost::array<double, 5>
{
public:
	pid_step_adjuster_coefficients()
	: boost::array<double, 5>({1.0/6, 2.0/6, 0, 0, 0})
	{}
};

template<>
class pid_step_adjuster_coefficients<adjuster_type::H0312> : public boost::array<double, 5>
{
public:
	pid_step_adjuster_coefficients()
	: boost::array<double, 5>({1.0/4, 2.0/2, 1.0/4, 3.0/4, 1.0/4})
	{}
};

template<>
class pid_step_adjuster_coefficients<adjuster_type::H312b> : public boost::array<double, 5>
{
public:
	pid_step_adjuster_coefficients()
	: boost::array<double, 5>({1.0/6, 2.0/6, 1.0/6, 3.0/6, 1.0/6})
	{}
};

template<>
class pid_step_adjuster_coefficients<adjuster_type::H312PID> : public boost::array<double, 5>
{
public:
	pid_step_adjuster_coefficients()
	: boost::array<double, 5>({1.0/18, 2.0/9, 1.0/18, 0, 0})
	{}
};

template<>
class pid_step_adjuster_coefficients<adjuster_type::H0321> : public boost::array<double, 5>
{
public:
	pid_step_adjuster_coefficients()
	: boost::array<double, 5>({5.0/4, 1.0/2, -3.0/4, -1.0/4, -3.0/4})
	{}
};

template<>
class pid_step_adjuster_coefficients<adjuster_type::H321> : public boost::array<double, 5>
{
public:
	pid_step_adjuster_coefficients()
	: boost::array<double, 5>({1.0/3, 1.0/18, -5.0/18, -5.0/6, -1.0/6})
	{}
};

}
}
}
}

#endif