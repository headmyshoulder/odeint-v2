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

#include <boost/numeric/odeint/algebra/standard_algebra.hpp>
#include <boost/numeric/odeint/algebra/standard_operations.hpp>
#include <boost/numeric/odeint/algebra/standard_resize.hpp>

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
    class Time = double ,
	class Algebra = standard_algebra ,
	class Operations = standard_operations< Time > ,
	class AdjustSizePolicy = adjust_size_initially_tag
	>
class explicit_error_rk54_ck
: public explicit_stepper_and_error_stepper_base<
	  explicit_error_rk54_ck< State , Time , Algebra , Operations , AdjustSizePolicy > ,
	  5 , 5 , 4 , State , Time , Algebra , Operations , AdjustSizePolicy >
{

public :

	BOOST_ODEINT_EXPLICIT_STEPPERS_AND_ERROR_STEPPERS_TYPEDEFS( explicit_error_rk54_ck , 5 , 5 , 4);

	typedef explicit_error_stepper_tag stepper_category;

	explicit_error_rk54_ck( void ) : m_size_adjuster() , m_x1() , m_x2() , m_x3() , m_x4() , m_x5() , m_x6()
	{
		boost::numeric::odeint::construct( m_x1 );
		boost::numeric::odeint::construct( m_x2 );
		boost::numeric::odeint::construct( m_x3 );
		boost::numeric::odeint::construct( m_x4 );
		boost::numeric::odeint::construct( m_x5 );
		boost::numeric::odeint::construct( m_x6 );
		m_size_adjuster.register_state( 0 , m_x1 );
		m_size_adjuster.register_state( 1 , m_x2 );
		m_size_adjuster.register_state( 2 , m_x3 );
		m_size_adjuster.register_state( 3 , m_x4 );
		m_size_adjuster.register_state( 4 , m_x5 );
		m_size_adjuster.register_state( 5 , m_x6 );
	}

	~explicit_error_rk54_ck( void )
	{
		boost::numeric::odeint::destruct( m_x1 );
		boost::numeric::odeint::destruct( m_x2 );
		boost::numeric::odeint::destruct( m_x3 );
		boost::numeric::odeint::destruct( m_x4 );
		boost::numeric::odeint::destruct( m_x5 );
		boost::numeric::odeint::destruct( m_x6 );
	}


	template< class System >
	void do_step_impl( System system , const state_type &in , const state_type &dxdt , const time_type t , state_type &out , const time_type dt , state_type &xerr )
	{
		const time_type c1 = static_cast<time_type> ( 37.0 ) / static_cast<time_type>( 378.0 );
		const time_type c3 = static_cast<time_type> ( 250.0 ) / static_cast<time_type>( 621.0 );
		const time_type c4 = static_cast<time_type> ( 125.0 ) / static_cast<time_type>( 594.0 );
		const time_type c6 = static_cast<time_type> ( 512.0 ) / static_cast<time_type>( 1771.0 );

		const time_type dc1 = c1 - static_cast<time_type> ( 2825.0 ) / static_cast<time_type>( 27648 );
		const time_type dc3 = c3 - static_cast<time_type> ( 18575.0 ) / static_cast<time_type>( 48384.0 );
		const time_type dc4 = c4 - static_cast<time_type> ( 13525.0 ) / static_cast<time_type>( 55296.0 );
		const time_type dc5 = static_cast<time_type> ( -277.0 ) / static_cast<time_type>( 14336.0 );
		const time_type dc6 = c6 - static_cast<time_type> ( 0.25 );

		do_step_impl( system , in , dxdt , t , out , dt );

		//error estimate
		algebra_type::for_each6( xerr , dxdt , m_x3 , m_x4 , m_x5 , m_x6 ,
					typename operations_type::scale_sum5( dt*dc1 , dt*dc3 , dt*dc4 , dt*dc5 , dt*dc6 ));

	}



	template< class System >
	void do_step_impl( System system , const state_type &in , const state_type &dxdt , const time_type t , state_type &out , const time_type dt )
	{
		const time_type a2 = static_cast<time_type> ( 0.2 );
		const time_type a3 = static_cast<time_type> ( 0.3 );
		const time_type a4 = static_cast<time_type> ( 0.6 );
		const time_type a5 = static_cast<time_type> ( 1.0 );
		const time_type a6 = static_cast<time_type> ( 0.875 );

		const time_type b21 = static_cast<time_type> ( 0.2 );
		const time_type b31 = static_cast<time_type> ( 3.0 ) / static_cast<time_type>( 40.0 );
		const time_type b32 = static_cast<time_type> ( 9.0 ) / static_cast<time_type>( 40.0 );
		const time_type b41 = static_cast<time_type> ( 0.3 );
		const time_type b42 = static_cast<time_type> ( -0.9 );
		const time_type b43 = static_cast<time_type> ( 1.2 );
		const time_type b51 = static_cast<time_type> ( -11.0 ) / static_cast<time_type>( 54.0 );
		const time_type b52 = static_cast<time_type> ( 2.5 );
		const time_type b53 = static_cast<time_type> ( -70.0 ) / static_cast<time_type>( 27.0 );
		const time_type b54 = static_cast<time_type> ( 35.0 ) / static_cast<time_type>( 27.0 );
		const time_type b61 = static_cast<time_type> ( 1631.0 ) / static_cast<time_type>( 55296.0 );
		const time_type b62 = static_cast<time_type> ( 175.0 ) / static_cast<time_type>( 512.0 );
		const time_type b63 = static_cast<time_type> ( 575.0 ) / static_cast<time_type>( 13824.0 );
		const time_type b64 = static_cast<time_type> ( 44275.0 ) / static_cast<time_type>( 110592.0 );
		const time_type b65 = static_cast<time_type> ( 253.0 ) / static_cast<time_type>( 4096.0 );

		const time_type c1 = static_cast<time_type> ( 37.0 ) / static_cast<time_type>( 378.0 );
		const time_type c3 = static_cast<time_type> ( 250.0 ) / static_cast<time_type>( 621.0 );
		const time_type c4 = static_cast<time_type> ( 125.0 ) / static_cast<time_type>( 594.0 );
		const time_type c6 = static_cast<time_type> ( 512.0 ) / static_cast<time_type>( 1771.0 );

		m_size_adjuster.adjust_size_by_policy( in , adjust_size_policy() );

		typename boost::unwrap_reference< System >::type &sys = system;

		//m_x1 = x + dt*b21*dxdt
		algebra_type::for_each3( m_x1 , in , dxdt ,
					typename operations_type::scale_sum2( 1.0 , dt*b21 ) );

		sys( m_x1 , m_x2 , t + dt*a2 );
		// m_x1 = x + dt*b31*dxdt + dt*b32*m_x2
		algebra_type::for_each4( m_x1 , in , dxdt , m_x2 ,
					typename operations_type::scale_sum3( 1.0 , dt*b31 , dt*b32 ));

		sys( m_x1 , m_x3 , t + dt*a3 );
		// m_x1 = x + dt * (b41*dxdt + b42*m_x2 + b43*m_x3)
		algebra_type::for_each5( m_x1 , in , dxdt , m_x2 , m_x3 ,
					typename operations_type::scale_sum4( 1.0 , dt*b41 , dt*b42 , dt*b43 ));

		sys( m_x1, m_x4 , t + dt*a4 );
		algebra_type::for_each6( m_x1 , in , dxdt , m_x2 , m_x3 , m_x4 ,
					typename operations_type::scale_sum5( 1.0 , dt*b51 , dt*b52 , dt*b53 , dt*b54 ));

		sys( m_x1 , m_x5 , t + dt*a5 );
		algebra_type::for_each7( m_x1 , in , dxdt , m_x2 , m_x3 , m_x4 , m_x5 ,
							typename operations_type::scale_sum6( 1.0 , dt*b61 , dt*b62 , dt*b63 , dt*b64 , dt*b65 ));

		sys( m_x1 , m_x6 , t + dt*a6 );
		algebra_type::for_each6( out , in , dxdt , m_x3 , m_x4 , m_x6 ,
					typename operations_type::scale_sum5( 1.0 , dt*c1 , dt*c3 , dt*c4 , dt*c6 ));

	}


	void adjust_size( const state_type &x )
	{
		m_size_adjuster.adjust_size( x );
		stepper_base_type::adjust_size( x );
	}


private:

    size_adjuster< state_type , 6 > m_size_adjuster;
    state_type m_x1, m_x2, m_x3, m_x4, m_x5, m_x6;

};







} // odeint
} // numeric
} // boost




#endif //BOOST_NUMERIC_ODEINT_STEPPER_EXPLICIT_ERROR_RK54_CK_HPP_INCLUDED
