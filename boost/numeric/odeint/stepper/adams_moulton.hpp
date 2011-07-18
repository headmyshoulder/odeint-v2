/*
 * adams_moulton.hpp
 *
 *  Created on: May 15, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_MOULTON_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_MOULTON_HPP_

#include <boost/ref.hpp>
#include <boost/bind.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>

#include <boost/numeric/odeint/stepper/detail/adams_moulton_call_algebra.hpp>
#include <boost/numeric/odeint/stepper/detail/adams_moulton_coefficients.hpp>
#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>






namespace boost {
namespace numeric {
namespace odeint {


/*
 * Static implicit Adams-Moulton multistep-solver without step size control and without dense output.
 *
 * # Define the number of steps
 */
template<
	size_t Steps ,
    class State ,
    class Value = double ,
    class Deriv = State ,
    class Time = Value ,
	class Algebra = range_algebra ,
	class Operations = default_operations ,
	class Resizer = initially_resizer
	>
class adams_moulton
{
private:


public :

	typedef State state_type;
	typedef state_wrapper< state_type > wrapped_state_type;
	typedef Value value_type;
	typedef Deriv deriv_type;
	typedef state_wrapper< deriv_type > wrapped_deriv_type;
	typedef Time time_type;
	typedef Algebra algebra_type;
	typedef Operations operations_type;
	typedef Resizer resizer_type;
	typedef stepper_tag stepper_category;

	typedef adams_moulton< Steps , State , Value , Deriv , Time , Algebra , Operations , Resizer > stepper_type;

	static const size_t steps = Steps;

	typedef unsigned short order_type;
	static const order_type order_value = steps + 1;

	typedef detail::rotating_buffer< wrapped_deriv_type , steps > step_storage_type;

    adams_moulton( ) : m_algebra( m_algebra_instance )
    { }

    adams_moulton( algebra_type &algebra ) : m_algebra( algebra )
    { }

	adams_moulton& operator=( const adams_moulton &stepper )
	{
		m_dxdt = stepper.m_dxdt;
		m_resizer = stepper.m_resizer;
		m_algebra = stepper.m_algebra;
		return *this;
	}

    order_type order( void ) const { return order_value; }


	/*
	 * Version 1 : do_step( system , x , t , dt , buf );
	 *
	 * solves the forwarding problem
	 */
    template< class System , class StateInOut , class ABBuf >
    void do_step( System system , StateInOut &in , const time_type &t , const time_type &dt , const ABBuf &buf )
    {
        do_step( system , in , t , in , dt , buf );
    }

    template< class System , class StateInOut , class ABBuf >
    void do_step( System system , const StateInOut &in , const time_type &t , const time_type &dt , const ABBuf &buf )
    {
        do_step( system , in , t , in , dt , buf );
    }



    /*
     * Version 2 : do_step( system , in , t , out , dt , buf );
     *
     * solves the forwarding problem
     */
    template< class System , class StateIn , class StateOut , class ABBuf >
    void do_step( System system , const StateIn &in , const time_type &t , StateOut &out , const time_type &dt , const ABBuf &buf )
    {
        typename boost::unwrap_reference< System >::type &sys = system;
        m_resizer.adjust_size( in , boost::bind( &stepper_type::resize<StateIn> , boost::ref( *this ) , _1 ) );
        sys( in , m_dxdt.m_v , t );
        detail::adams_moulton_call_algebra< steps , algebra_type , operations_type >()( m_algebra , in , out , m_dxdt.m_v , buf , m_coefficients , dt );
    }

	template< class System , class StateIn , class StateOut , class ABBuf >
	void do_step( System system , const StateIn &in , const time_type &t , const StateOut &out , const time_type &dt , const ABBuf &buf )
	{
		typename boost::unwrap_reference< System >::type &sys = system;
        m_resizer.adjust_size( in , boost::bind( &stepper_type::resize<StateIn> , boost::ref( *this ) , _1 ) );
		sys( in , m_dxdt.m_v , t );
		detail::adams_moulton_call_algebra< steps , algebra_type , operations_type >()( m_algebra , in , out , m_dxdt.m_v , buf , m_coefficients , dt );
	}



	template< class StateIn >
	bool resize( const StateIn &x )
	{
	    return adjust_size_by_resizeability( m_dxdt , x , typename wrapped_deriv_type::is_resizeable() );
	}


	template< class StateType >
	void adjust_size( const StateType &x )
	{
		resize( x );
	}

    algebra_type& algebra()
    {   return m_algebra; }

    const algebra_type& algebra() const
    {   return m_algebra; }


private:

	const detail::adams_moulton_coefficients< value_type , steps > m_coefficients;
	wrapped_deriv_type m_dxdt;
	resizer_type m_resizer;

protected:
	algebra_type m_algebra_instance;
	algebra_type &m_algebra;
};




} // odeint
} // numeric
} // boost



#endif /* BOOST_NUMERIC_ODEINT_STEPPER_ADAMS_MOULTON_HPP_ */
