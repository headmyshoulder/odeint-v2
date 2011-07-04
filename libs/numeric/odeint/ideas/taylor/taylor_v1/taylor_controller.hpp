/*
 boost header: NUMERIC_ODEINT_STEPPER/controlled_error_stepper.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_TAYLOR_CONTROLLER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_TAYLOR_CONTROLLER_HPP_INCLUDED

#include <cmath>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>




namespace boost {
namespace numeric {
namespace odeint {


/*
 * explicit stepper version
 */
template
<
	class TaylorStepper
>
class taylor_controller
{

public:

	typedef TaylorStepper stepper_type;
	typedef typename stepper_type::state_type state_type;
	typedef typename stepper_type::value_type value_type;
	typedef typename stepper_type::deriv_type deriv_type;
	typedef typename stepper_type::time_type time_type;
	typedef typename stepper_type::order_type order_type;
	typedef default_error_checker< value_type > error_checker_type;
	typedef controlled_stepper_tag stepper_category;


	taylor_controller(
			value_type eps_abs = static_cast< value_type >( 1.0e-6 ) ,
			value_type eps_rel = static_cast< value_type >( 1.0e-6 ) ,
			value_type a_x = static_cast< value_type >( 1.0 ) ,
			value_type a_dxdt = static_cast< value_type >( 1.0 ) ,
			const stepper_type &stepper = stepper_type()
			)
	: m_stepper( stepper ) ,
	  m_error_checker( eps_abs , eps_rel , a_x , a_dxdt ) ,
	  m_xerr() , m_xnew()
	{
	}


	template< class System >
	controlled_step_result try_step( System system , state_type &x , time_type &t , time_type &dt )
	{
		controlled_step_result res = try_step( system , x , t , m_xnew , dt );
		if( ( res == success_step_size_increased ) || ( res == success_step_size_unchanged ) )
		{
			x = m_xnew;
		}
		return res;
	}


	template< class System >
	controlled_step_result try_step( System system , state_type &in , time_type &t , state_type &out , time_type &dt )
	{
		using std::max;
		using std::min;
		using std::pow;

		// do one step with error calculation
		m_stepper.do_step( system , in , t , out , dt , m_xerr );

		value_type max_rel_error = m_error_checker.error( in , m_stepper.get_last_derivs()[0] , m_xerr , dt );

		if( max_rel_error > 1.1 )
		{
			// error too large - decrease dt ,limit scaling factor to 0.2 and reset state
			dt *= max( 0.9 * pow( max_rel_error , -1.0 / ( m_stepper.error_order() - 1.0 ) ) , 0.2 );
			return step_size_decreased;
		}
		else
		{
			if( max_rel_error < 0.5 )
			{
				//error too small - increase dt and keep the evolution and limit scaling factor to 5.0
				t += dt;
				dt *= min( 0.9 * pow( max_rel_error , -1.0 / m_stepper.stepper_order() ) , 5.0 );
				return success_step_size_increased;
			}
			else
			{
				t += dt;
				return success_step_size_unchanged;
			}
		}
	}


	stepper_type& stepper( void )
	{
		return m_stepper;
	}

	const stepper_type& stepper( void ) const
	{
		return m_stepper;
	}


private:



	stepper_type m_stepper;
	error_checker_type m_error_checker;
	state_type m_xerr;
	state_type m_xnew;
};








} // odeint
} // numeric
} // boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_TAYLOR_CONTROLLER_HPP_INCLUDED
