/*
 boost header: NUMERIC_ODEINT_STEPPER/explicit_rk4.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_EXPLICIT_RK4_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_EXPLICIT_RK4_HPP_INCLUDED

#include <boost/ref.hpp>

#include <boost/numeric/odeint/algebra/standard_algebra.hpp>
#include <boost/numeric/odeint/algebra/standard_operations.hpp>
#include <boost/numeric/odeint/algebra/standard_resize.hpp>

#include <boost/numeric/odeint/stepper/base/explicit_stepper_base.hpp>
#include <boost/numeric/odeint/stepper/detail/macros.hpp>

namespace boost {
namespace numeric {
namespace odeint {


template<
    class State ,
    class Value = double ,
    class Deriv = State ,
    class Time = Value ,
	class Algebra = standard_algebra ,
	class Operations = standard_operations ,
	class AdjustSizePolicy = adjust_size_initially_tag
	>
class explicit_rk4
: public explicit_stepper_base<
	  explicit_rk4< State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > ,
	  4 , State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy >
{
	void initialize( void )
	{
		boost::numeric::odeint::construct( m_dxt );
		boost::numeric::odeint::construct( m_dxm );
		boost::numeric::odeint::construct( m_dxh );
		boost::numeric::odeint::construct( m_x_tmp );
		m_deriv_adjuster.register_state( 0 , m_dxt );
		m_deriv_adjuster.register_state( 1 , m_dxm );
		m_deriv_adjuster.register_state( 2 , m_dxh );
		m_state_adjuster.register_state( 0 , m_x_tmp );
	}

	void copy( const explicit_rk4 &rk )
	{
		boost::numeric::odeint::copy( rk.m_dxt , m_dxt );
		boost::numeric::odeint::copy( rk.m_dxm , m_dxm );
		boost::numeric::odeint::copy( rk.m_dxh , m_dxh );
		boost::numeric::odeint::copy( rk.m_x_tmp , m_x_tmp );
	}

public :


	BOOST_ODEINT_EXPLICIT_STEPPERS_TYPEDEFS( explicit_rk4 , 4 );


	explicit_rk4( void )
	: stepper_base_type() , m_deriv_adjuster() , m_state_adjuster() , m_dxt() , m_dxm() , m_dxh() , m_x_tmp()
	{
		initialize();
	}

	~explicit_rk4( void )
	{
		boost::numeric::odeint::destruct( m_dxt );
		boost::numeric::odeint::destruct( m_dxm );
		boost::numeric::odeint::destruct( m_dxh );
		boost::numeric::odeint::destruct( m_x_tmp );
	}

	explicit_rk4( const explicit_rk4 &rk )
	: stepper_base_type( rk ) , m_deriv_adjuster() , m_state_adjuster() , m_dxt() , m_dxm() , m_dxh() , m_x_tmp()
	{
		initialize();
		copy( rk );
	}

	explicit_rk4& operator=( const explicit_rk4 &rk )
	{
		stepper_base_type::operator=( rk );
		copy( rk );
		return *this;
	}


	template< class System , class StateIn , class DerivIn , class StateOut >
	void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt , const time_type &t , StateOut &out , const time_type &dt )
	{
		// ToDo : check if size of in,dxdt,out are equal?

        const value_type val1 = static_cast< value_type >( 1.0 );

		m_deriv_adjuster.adjust_size_by_policy( in , adjust_size_policy() );
		m_state_adjuster.adjust_size_by_policy( in , adjust_size_policy() );

		typename boost::unwrap_reference< System >::type &sys = system;

        time_type dh = static_cast< value_type >( 0.5 ) * dt;
        time_type th = t + dh;

        // dt * dxdt = k1
        // m_x_tmp = x + dh*dxdt
        typename algebra_type::for_each3()( m_x_tmp , in , dxdt ,
        		typename operations_type::template scale_sum2< value_type , time_type >( val1 , dh ) );


        // dt * m_dxt = k2
        sys( m_x_tmp , m_dxt , th );

        // m_x_tmp = x + dh*m_dxt
        typename algebra_type::for_each3()( m_x_tmp , in , m_dxt ,
        		typename operations_type::template scale_sum2< value_type , time_type >( val1 , dh ) );


        // dt * m_dxm = k3
        sys( m_x_tmp , m_dxm , th );
        //m_x_tmp = x + dt*m_dxm
        typename algebra_type::for_each3()( m_x_tmp , in , m_dxm ,
        		typename operations_type::template scale_sum2< value_type , time_type >( val1 , dt ) );


        // dt * m_dxh = k4
        sys( m_x_tmp , m_dxh , t + dt );
        //x += dt/6 * ( m_dxdt + m_dxt + val2*m_dxm )
        time_type dt6 = dt / static_cast< value_type >( 6.0 );
        time_type dt3 = dt / static_cast< value_type >( 3.0 );
        typename algebra_type::for_each6()( out , in , dxdt , m_dxt , m_dxm , m_dxh ,
        		typename operations_type::template scale_sum5< value_type , time_type , time_type , time_type , time_type >( 1.0 , dt6 , dt3 , dt3 , dt6 ) );
	}


	template< class StateType >
	void adjust_size( const StateType &x )
	{
		m_deriv_adjuster.adjust_size( x );
		m_state_adjuster.adjust_size( x );
		stepper_base_type::adjust_size( x );
	}


private:

	size_adjuster< deriv_type , 3 > m_deriv_adjuster;
	size_adjuster< state_type , 1 > m_state_adjuster;

    deriv_type m_dxt;
    deriv_type m_dxm;
    deriv_type m_dxh;
    state_type m_x_tmp;

};




} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_STEPPER_EXPLICIT_RK4_HPP_INCLUDED
