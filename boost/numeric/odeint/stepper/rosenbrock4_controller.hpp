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


namespace boost {
namespace numeric {
namespace odeint {

template< class Stepper >
class rosenbrock4_controller
{
public:

	typedef Stepper stepper_type;
	typedef typename stepper_type::value_type value_type;
	typedef typename stepper_type::state_type state_type;
	typedef typename stepper_type::time_type time_type;
	typedef typename stepper_type::deriv_type deriv_type;
	typedef controlled_stepper_tag stepper_category;


	rosenbrock4_controller( value_type atol = 1.0e-6 , value_type rtol = 1.0e-6 , const stepper_type &stepper = stepper_type() )
    : m_stepper() , m_atol( atol ) , m_rtol( rtol ) ,
      m_first_step( true ) , m_err_old( 0.0 ) , m_dt_old( 0.0 ) ,
      m_last_rejected( false )
	{
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
	 * ToDo : xerr in size_adjuster
	 */
	template< class System >
	boost::numeric::odeint::controlled_step_result
	try_step( System sys , const state_type &x , value_type &t , state_type &xout , value_type &dt )
	{
		static const value_type safe = 0.9 , fac1 = 5.0 , fac2 = 1.0 / 6.0;

		const size_t n = x.size();
		state_type xerr( n );
		m_stepper.do_step( sys , x , t , xout , dt , xerr );
		value_type err = error( xout , x , xerr );

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
	value_type m_atol , m_rtol;
	bool m_first_step;
	value_type m_err_old , m_dt_old;
	bool m_last_rejected;
};






} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_STEPPER_ROSENBROCK4_CONTROLLER_HPP_ */
