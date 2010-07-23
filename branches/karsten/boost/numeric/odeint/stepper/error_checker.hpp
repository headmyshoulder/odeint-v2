/*
 boost header: NUMERIC_ODEINT/error_checker.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_ERROR_CHECKER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_ERROR_CHECKER_HPP_INCLUDED

namespace boost {
namespace numeric {
namespace odeint {


template< class State , class Time >
class error_checker_standard
{
public:

	typedef State state_type;
	typedef Time time_type;


	error_checker_standard( void )
	{}

	time_type error( const state_type &x_old , const state_type &dxdt_old , const state_type &x_err , time_type dt )
	{
		return 0.0;
	}

private:

	//	time_type m_eps_abs;
	//	time_type m_eps_rel;
	//	time_type m_a_x;
	//	time_type m_a_dxdt;
};

} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_ERROR_CHECKER_HPP_INCLUDED
