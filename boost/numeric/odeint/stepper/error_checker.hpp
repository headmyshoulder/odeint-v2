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

#include <boost/numeric/odeint/algebra/standard_algebra.hpp>
#include <boost/numeric/odeint/algebra/standard_operations.hpp>

namespace boost {
namespace numeric {
namespace odeint {


template
<
	class State ,
	class Time ,
	class Algebra = standard_algebra ,
	class Operations = standard_operations
>
class error_checker_standard
{
public:

	typedef State state_type;
	typedef Time time_type;
	typedef Algebra algebra_type;
	typedef Operations operations_type;


	error_checker_standard( void ) : m_eps_abs( 1E-6 ) , m_eps_rel( 1E-6 ) , m_a_x( 1.0 ) , m_a_dxdt( 1.0 )
	{}


	/*
	 * ToDo: implement constructor with epsilons
	 */
	time_type error( const state_type &x_old , const state_type &dxdt_old , state_type &x_err , const time_type &dt )
	{
		// this overwrites x_err !
		algebra_type::for_each3( x_old , dxdt_old , x_err ,
					             typename operations_type::template rel_error< time_type >( m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt*dt ) );

		return algebra_type::template reduce< time_type >( x_err , typename operations_type::template maximum< time_type >() , 0.0 );
	}

private:

	time_type m_eps_abs;
	time_type m_eps_rel;
	time_type m_a_x;
	time_type m_a_dxdt;
};

} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_ERROR_CHECKER_HPP_INCLUDED
