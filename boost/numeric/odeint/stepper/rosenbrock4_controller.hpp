/*
 * rosenbrock4_controller.hpp
 *
 *  Created on: Jan 31, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_ROSENBROCK4_CONTROLLER_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_ROSENBROCK4_CONTROLLER_HPP_

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

#include <boost/numeric/odeint/stepper/rosenbrock4.hpp>

#include <iostream>


namespace boost {
namespace numeric {
namespace odeint {

template< class Stepper >
class rosenbrock4_controller
{
private:

	void initialize_variables( void )
	{
		m_state_adjuster.register_state( 0 , m_xerr );
	}

	void copy_variables( const rosenbrock4_controller &rb )
	{
        /** @ToDo: MSVC 10 fails on these copy operations with invalid null pointer exception, find out why
         * maybe add rosenbrock4 to stepper_copying test cases
         */
		m_stepper = rb.m_stepper;
		m_xerr = rb.m_xerr;
        /* end MSVC 10 error */
		m_atol = rb.m_atol;
		m_rtol = rb.m_rtol;
		m_first_step = rb.m_first_step;
		m_err_old = rb.m_err_old;
		m_dt_old = rb.m_dt_old;
		m_last_rejected = rb.m_last_rejected;
	}

public:

	typedef Stepper stepper_type;
	typedef typename stepper_type::value_type value_type;
	typedef typename stepper_type::state_type state_type;
	typedef typename stepper_type::time_type time_type;
	typedef typename stepper_type::deriv_type deriv_type;
	typedef typename stepper_type::adjust_size_policy adjust_size_policy;
	typedef controlled_stepper_tag stepper_category;


	rosenbrock4_controller( value_type atol = 1.0e-6 , value_type rtol = 1.0e-6 , const stepper_type &stepper = stepper_type() )
    : m_stepper() , m_state_adjuster() , m_xerr() ,
      m_atol( atol ) , m_rtol( rtol ) ,
      m_first_step( true ) , m_err_old( 0.0 ) , m_dt_old( 0.0 ) ,
      m_last_rejected( false )
	{
		initialize_variables();
	}

	rosenbrock4_controller( const rosenbrock4_controller &rb )
	: m_stepper() , m_state_adjuster() , m_xerr() ,
      m_atol( 1.0e-6l ) , m_rtol( 1.0e-6 ) ,
      m_first_step( true ) , m_err_old( 0.0 ) , m_dt_old( 0.0 ) ,
      m_last_rejected( false )
	{
        initialize_variables();
		copy_variables( rb );
	}

	rosenbrock4_controller& operator=( const rosenbrock4_controller &rb )
	{
		copy_variables( rb );
		return *this;
	}

	value_type error( const state_type &x , const state_type &xold , const state_type &xerr )
	{
		const size_t n = x.size();
		value_type err = 0.0 , sk = 0.0;
		for( size_t i=0 ; i<n ; ++i )
		{
			sk = m_atol + m_rtol * std::max( std::abs( xold[i] ) , std::abs( x[i] ) );
			err += xerr[i] * xerr[i] / sk / sk;
		}
		return sqrt( err / value_type( n ) );
	}

	value_type last_error( void ) const
	{
		return m_err_old;
	}




	/*
	 * ToDo : xout in size_adjuster
	 * Check
	 */
	template< class System >
	boost::numeric::odeint::controlled_step_result
	try_step( System sys , state_type &x , value_type &t , value_type &dt )
	{
		state_type xout( x.size() );
		boost::numeric::odeint::controlled_step_result res = try_step( sys , x , t , xout , dt );
		x = xout;
		return res;
	}


	/*
	 * ToDo : xerr in size_adjuster und adjust size
	 */
	template< class System >
	boost::numeric::odeint::controlled_step_result
	try_step( System sys , const state_type &x , value_type &t , state_type &xout , value_type &dt )
	{
		static const value_type safe = 0.9 , fac1 = 5.0 , fac2 = 1.0 / 6.0;

		m_state_adjuster.adjust_size_by_policy( x , adjust_size_policy() );

		m_stepper.do_step( sys , x , t , xout , dt , m_xerr );
		value_type err = error( xout , x , m_xerr );

		value_type fac = std::max( fac2 ,std::min( fac1 , std::pow( err , 0.25 ) / safe ) );
		value_type dt_new = dt / fac;
		if ( err <= 1.0 )
		{
			if( m_first_step )
			{
				m_first_step = false;
			}
			else
			{
				value_type fac_pred = ( m_dt_old / dt ) * pow( err * err / m_err_old , 0.25 ) / safe;
				fac_pred = std::max( fac2 , std::min( fac1 , fac_pred ) );
				fac = std::max( fac , fac_pred );
				dt_new = dt / fac;
			}

			m_dt_old = dt;
			m_err_old = std::max( 0.01 , err );
			if( m_last_rejected )
				dt_new = ( dt >= 0.0 ? std::min( dt_new , dt ) : std::max( dt_new , dt ) );
			t += dt;
			dt = dt_new;
			m_last_rejected = false;
			return success_step_size_increased;
		}
		else
		{
			dt = dt_new;
			m_last_rejected = true;
			return step_size_decreased;
		}
	}


	template< class StateType >
	void adjust_size( const StateType &x )
	{
		m_stepper.adjust_size( x );
		m_state_adjuster.adjust_size( x );
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
	size_adjuster< state_type , 1 > m_state_adjuster;
	state_type m_xerr;
	value_type m_atol , m_rtol;
	bool m_first_step;
	value_type m_err_old , m_dt_old;
	bool m_last_rejected;
};






} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_STEPPER_ROSENBROCK4_CONTROLLER_HPP_ */
