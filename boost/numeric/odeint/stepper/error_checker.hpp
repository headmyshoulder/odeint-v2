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


/*
 * ToDo: implement constructor with epsilons
 */
template
<
	class Value ,
	class Algebra = standard_algebra ,
	class Operations = standard_operations
>
class error_checker_standard
{
public:

	typedef Value value_type;
	typedef Algebra algebra_type;
	typedef Operations operations_type;


	error_checker_standard( void ) : m_eps_abs( 1E-6 ) , m_eps_rel( 1E-6 ) , m_a_x( 1.0 ) , m_a_dxdt( 1.0 )
	{}


	template< class State , class Deriv , class Err , class Time >
	value_type error( const State &x_old , const Deriv &dxdt_old , Err &x_err , const Time &dt )
	{
		// this overwrites x_err !
		typename algebra_type::for_each3()( x_old , dxdt_old , x_err ,
					             typename operations_type::template rel_error< value_type >( m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt * detail::get_value( dt ) ) );

		value_type res = typename algebra_type::reduce()( x_err , typename operations_type::template maximum< value_type >() , 0.0 );
		return res;
	}

private:

	value_type m_eps_abs;
	value_type m_eps_rel;
	value_type m_a_x;
	value_type m_a_dxdt;
};

} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_ERROR_CHECKER_HPP_INCLUDED
