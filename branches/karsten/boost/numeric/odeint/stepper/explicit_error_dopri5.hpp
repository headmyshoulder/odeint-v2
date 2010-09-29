/*
 boost header: boost_numeric_odeint/explicit_error_dopri5.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_STEPPER_EXPLICIT_ERROR_DOPRI5_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_STEPPER_EXPLICIT_ERROR_DOPRI5_HPP_INCLUDED

#include <boost/numeric/odeint/algebra/standard_algebra.hpp>
#include <boost/numeric/odeint/algebra/standard_operations.hpp>
#include <boost/numeric/odeint/algebra/standard_resize.hpp>

#include <boost/numeric/odeint/stepper/base/explicit_error_stepper_base.hpp>
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
	class Algebra = standard_algebra< State > ,
	class Operations = standard_operations< Time > ,
	class AdjustSizePolicy = adjust_size_initially_tag
	>
class explicit_error_dopri5
: public explicit_error_stepper_base<
	  explicit_error_dopri5< State , Time , Algebra , Operations , AdjustSizePolicy > ,
	  5 , 4 , State , Time , Algebra , Operations , AdjustSizePolicy >
{

public :

	BOOST_ODEINT_EXPLICIT_ERROR_STEPPERS_TYPEDEFS( explicit_error_dopri5 , 5 , 4 );

	explicit_error_dopri5( void )
	: m_size_adjuster() , m_x1() , m_x2() , m_x3() , m_x4() , m_x5() , m_x6()
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

	~explicit_error_dopri5( void )
	{
		boost::numeric::odeint::destruct( m_x1 );
		boost::numeric::odeint::destruct( m_x2 );
		boost::numeric::odeint::destruct( m_x3 );
		boost::numeric::odeint::destruct( m_x4 );
		boost::numeric::odeint::destruct( m_x5 );
		boost::numeric::odeint::destruct( m_x6 );
	}


	template< class System >
	void do_step_impl( System system , const state_type &in , const state_type &dxdt , state_type &out , time_type t , time_type dt , state_type &xerr )
	{
	}



	template< class System >
	void do_step_impl( System system , const state_type &in , const state_type &dxdt , state_type &out , time_type t , time_type dt )
	{
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

#endif //BOOST_BOOST_NUMERIC_ODEINT_STEPPER_EXPLICIT_ERROR_DOPRI5_HPP_INCLUDED
