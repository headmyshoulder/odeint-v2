/*
 boost header: NUMERIC_ODEINT_STEPPER/explicit_error_rk54_ck.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_EXPLICIT_ERROR_RK54_CK_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_EXPLICIT_ERROR_RK54_CK_HPP_INCLUDED

#include <boost/ref.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/algebra/default_resize.hpp>

#include <boost/numeric/odeint/stepper/base/explicit_stepper_and_error_stepper_base.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/detail/macros.hpp>

namespace boost {
namespace numeric {
namespace odeint {





/*
 * ToDo: Check orders rk_ckc
 */
template<
	class State ,
    class Value = double ,
    class Deriv = State ,
    class Time = double ,
	class Algebra = range_algebra ,
	class Operations = default_operations ,
	class AdjustSizePolicy = adjust_size_initially_tag
	>
class explicit_error_rk54_ck
: public explicit_stepper_and_error_stepper_base<
	  explicit_error_rk54_ck< State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > ,
	  5 , 5 , 4 , State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy >
{

	void initialize( void )
	{
		boost::numeric::odeint::construct( m_x_tmp );
		boost::numeric::odeint::construct( m_k2 );
		boost::numeric::odeint::construct( m_k3 );
		boost::numeric::odeint::construct( m_k4 );
		boost::numeric::odeint::construct( m_k5 );
		boost::numeric::odeint::construct( m_k6 );
		m_state_adjuster.register_state( 0 , m_x_tmp );
		m_deriv_adjuster.register_state( 0 , m_k2 );
		m_deriv_adjuster.register_state( 1 , m_k3 );
		m_deriv_adjuster.register_state( 2 , m_k4 );
		m_deriv_adjuster.register_state( 3 , m_k5 );
		m_deriv_adjuster.register_state( 4 , m_k6 );
	}

	void copy( const explicit_error_rk54_ck &rk )
	{
		boost::numeric::odeint::copy( rk.m_x_tmp , m_x_tmp );
		boost::numeric::odeint::copy( rk.m_k2 , m_k2 );
		boost::numeric::odeint::copy( rk.m_k3 , m_k3 );
		boost::numeric::odeint::copy( rk.m_k4 , m_k4 );
		boost::numeric::odeint::copy( rk.m_k5 , m_k5 );
		boost::numeric::odeint::copy( rk.m_k6 , m_k6 );
	}

public :

	BOOST_ODEINT_EXPLICIT_STEPPERS_AND_ERROR_STEPPERS_TYPEDEFS( explicit_error_rk54_ck , 5 , 5 , 4);

	typedef explicit_error_stepper_tag stepper_category;

	explicit_error_rk54_ck( void )
	: stepper_base_type() , m_state_adjuster() , m_deriv_adjuster() , m_x_tmp() , m_k2() , m_k3() , m_k4() , m_k5() , m_k6()
	{
		initialize();
	}

	~explicit_error_rk54_ck( void )
	{
		boost::numeric::odeint::destruct( m_x_tmp );
		boost::numeric::odeint::destruct( m_k2 );
		boost::numeric::odeint::destruct( m_k3 );
		boost::numeric::odeint::destruct( m_k4 );
		boost::numeric::odeint::destruct( m_k5 );
		boost::numeric::odeint::destruct( m_k6 );
	}

	explicit_error_rk54_ck( const explicit_error_rk54_ck &rk )
	: stepper_base_type( rk ) , m_state_adjuster() , m_deriv_adjuster() , m_x_tmp() , m_k2() , m_k3() , m_k4() , m_k5() , m_k6()
	{
		initialize();
		copy( rk );
	}

	explicit_error_rk54_ck& operator=( const explicit_error_rk54_ck &rk )
	{
		stepper_base_type::operator=( rk );
		copy( rk );
		return *this;
	}

	template< class System , class StateIn , class DerivIn , class StateOut , class Err >
	void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt , const time_type &t , StateOut &out , const time_type &dt , Err &xerr )
	{
		const value_type c1 = static_cast<value_type> ( 37.0 ) / static_cast<value_type>( 378.0 );
		const value_type c3 = static_cast<value_type> ( 250.0 ) / static_cast<value_type>( 621.0 );
		const value_type c4 = static_cast<value_type> ( 125.0 ) / static_cast<value_type>( 594.0 );
		const value_type c6 = static_cast<value_type> ( 512.0 ) / static_cast<value_type>( 1771.0 );

		const value_type dc1 = c1 - static_cast<value_type> ( 2825.0 ) / static_cast<value_type>( 27648 );
		const value_type dc3 = c3 - static_cast<value_type> ( 18575.0 ) / static_cast<value_type>( 48384.0 );
		const value_type dc4 = c4 - static_cast<value_type> ( 13525.0 ) / static_cast<value_type>( 55296.0 );
		const value_type dc5 = static_cast<value_type> ( -277.0 ) / static_cast<value_type>( 14336.0 );
		const value_type dc6 = c6 - static_cast<value_type> ( 0.25 );

		do_step_impl( system , in , dxdt , t , out , dt );

		//error estimate
		typename algebra_type::for_each6()( xerr , dxdt , m_k3 , m_k4 , m_k5 , m_k6 ,
				typename operations_type::template scale_sum5< time_type , time_type , time_type , time_type , time_type >( dt*dc1 , dt*dc3 , dt*dc4 , dt*dc5 , dt*dc6 ));

	}



	template< class System , class StateIn , class DerivIn , class StateOut >
	void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt , const time_type &t , StateOut &out , const time_type &dt )
	{
		const value_type a2 = static_cast<value_type> ( 0.2 );
		const value_type a3 = static_cast<value_type> ( 0.3 );
		const value_type a4 = static_cast<value_type> ( 0.6 );
		const value_type a5 = static_cast<value_type> ( 1.0 );
		const value_type a6 = static_cast<value_type> ( 0.875 );

		const value_type b21 = static_cast<value_type> ( 0.2 );
		const value_type b31 = static_cast<value_type> ( 3.0 ) / static_cast<value_type>( 40.0 );
		const value_type b32 = static_cast<value_type> ( 9.0 ) / static_cast<value_type>( 40.0 );
		const value_type b41 = static_cast<value_type> ( 0.3 );
		const value_type b42 = static_cast<value_type> ( -0.9 );
		const value_type b43 = static_cast<value_type> ( 1.2 );
		const value_type b51 = static_cast<value_type> ( -11.0 ) / static_cast<value_type>( 54.0 );
		const value_type b52 = static_cast<value_type> ( 2.5 );
		const value_type b53 = static_cast<value_type> ( -70.0 ) / static_cast<value_type>( 27.0 );
		const value_type b54 = static_cast<value_type> ( 35.0 ) / static_cast<value_type>( 27.0 );
		const value_type b61 = static_cast<value_type> ( 1631.0 ) / static_cast<value_type>( 55296.0 );
		const value_type b62 = static_cast<value_type> ( 175.0 ) / static_cast<value_type>( 512.0 );
		const value_type b63 = static_cast<value_type> ( 575.0 ) / static_cast<value_type>( 13824.0 );
		const value_type b64 = static_cast<value_type> ( 44275.0 ) / static_cast<value_type>( 110592.0 );
		const value_type b65 = static_cast<value_type> ( 253.0 ) / static_cast<value_type>( 4096.0 );

		const value_type c1 = static_cast<value_type> ( 37.0 ) / static_cast<value_type>( 378.0 );
		const value_type c3 = static_cast<value_type> ( 250.0 ) / static_cast<value_type>( 621.0 );
		const value_type c4 = static_cast<value_type> ( 125.0 ) / static_cast<value_type>( 594.0 );
		const value_type c6 = static_cast<value_type> ( 512.0 ) / static_cast<value_type>( 1771.0 );

		m_state_adjuster.adjust_size_by_policy( in , adjust_size_policy() );
		m_deriv_adjuster.adjust_size_by_policy( in , adjust_size_policy() );

		typename boost::unwrap_reference< System >::type &sys = system;

		//m_x1 = x + dt*b21*dxdt
		typename algebra_type::for_each3()( m_x_tmp , in , dxdt ,
				typename operations_type::template scale_sum2< value_type , time_type >( 1.0 , dt*b21 ) );

		sys( m_x_tmp , m_k2 , t + dt*a2 );
		// m_x_tmp = x + dt*b31*dxdt + dt*b32*m_x2
		typename algebra_type::for_each4()( m_x_tmp , in , dxdt , m_k2 ,
				typename operations_type::template scale_sum3< value_type , time_type , time_type >( 1.0 , dt*b31 , dt*b32 ));

		sys( m_x_tmp , m_k3 , t + dt*a3 );
		// m_x_tmp = x + dt * (b41*dxdt + b42*m_x2 + b43*m_x3)
		typename algebra_type::for_each5()( m_x_tmp , in , dxdt , m_k2 , m_k3 ,
				typename operations_type::template scale_sum4< value_type , time_type , time_type , time_type >( 1.0 , dt*b41 , dt*b42 , dt*b43 ));

		sys( m_x_tmp, m_k4 , t + dt*a4 );
		typename algebra_type::for_each6()( m_x_tmp , in , dxdt , m_k2 , m_k3 , m_k4 ,
				typename operations_type::template scale_sum5< value_type , time_type , time_type , time_type , time_type >( 1.0 , dt*b51 , dt*b52 , dt*b53 , dt*b54 ));

		sys( m_x_tmp , m_k5 , t + dt*a5 );
		typename algebra_type::for_each7()( m_x_tmp , in , dxdt , m_k2 , m_k3 , m_k4 , m_k5 ,
				typename operations_type::template scale_sum6< value_type , time_type , time_type , time_type , time_type , time_type >( 1.0 , dt*b61 , dt*b62 , dt*b63 , dt*b64 , dt*b65 ));

		sys( m_x_tmp , m_k6 , t + dt*a6 );
		typename algebra_type::for_each6()( out , in , dxdt , m_k3 , m_k4 , m_k6 ,
				typename operations_type::template scale_sum5< value_type , time_type , time_type , time_type , time_type >( 1.0 , dt*c1 , dt*c3 , dt*c4 , dt*c6 ));

	}


	template< class StateType >
	void adjust_size( const StateType &x )
	{
		m_state_adjuster.adjust_size( x );
		m_deriv_adjuster.adjust_size( x );
		stepper_base_type::adjust_size( x );
	}


private:

	size_adjuster< state_type , 1 > m_state_adjuster;
    size_adjuster< deriv_type , 5 > m_deriv_adjuster;
    state_type m_x_tmp;
    deriv_type m_k2, m_k3, m_k4, m_k5, m_k6;

};







} // odeint
} // numeric
} // boost




#endif //BOOST_NUMERIC_ODEINT_STEPPER_EXPLICIT_ERROR_RK54_CK_HPP_INCLUDED
