/*
 * controlled_error_bs.hpp
 *
 *  Created on: Jul 12, 2011
 *      Author: mario
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_ERROR_BS_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_ERROR_BS_HPP_

#include <boost/ref.hpp>

#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/stepper/detail/macros.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template<
    class State ,
    class Value = double ,
    class Deriv = State ,
    class Time = Value ,
    class Algebra = range_algebra ,
    class Operations = default_operations ,
    class Resizer = initially_resizer
    >
class controlled_error_bs {

public:

    typedef State state_type;
    typedef Value value_type;
    typedef Deriv deriv_type;
    typedef Time time_type;
    typedef Algebra algebra_type;
    typedef Operations operations_type;
    typedef Resizer resizer_type;
    typedef state_wrapper< state_type > wrapped_state_type;
    typedef state_wrapper< deriv_type > wrapped_deriv_type;

    typedef controlled_error_bs< State , Value , Deriv , Time , Algebra , Operations , Resizer > controlled_error_bs_type;


    controlled_error_bs(
            time_type abs_err = 1E-6 , time_type rel_err = 1E-6 ,
            time_type factor_x = 1.0 , time_type factor_dxdt = 1.0 )
        : m_error_checker( abs_err, rel_err, factor_x, factor_dxdt ),
          m_k_max(8) ,
          m_safety1(0.25) , m_safety2(0.7),
          m_max_dt_factor( 0.1 ) , m_min_step_scale(5E-5) , m_max_step_scale(0.7),
          m_continuous_calls(false) , m_decreased_step_during_last_call( false ) ,
          m_dt_last( 1.0E30) ,
          m_current_eps( -1.0 ) ,
          m_error( m_k_max ) ,
          m_a( m_k_max+1 ) ,
          m_alpha( m_k_max , value_vector( m_k_max ) ) ,
          m_interval_sequence( m_k_max+1 ) ,
          m_times( m_k_max ) ,
          m_d( m_k_max )
    {
        for( unsigned short i = 1; i <= m_k_max+1; i++ )
            m_interval_sequence[i-1] = (2*i);
    }

    /*
	 * Version 1 : try_step( sys , x , t , dt )
	 *
	 * The overloads are needed to solve the forwarding problem
	 */
	template< class System , class StateInOut >
	controlled_step_result try_step( System system , StateInOut &x , time_type &t , time_type &dt )
	{
		return try_step_v1( system , x , t, dt );
	}

	template< class System , class StateInOut >
	controlled_step_result try_step( System system , const StateInOut &x , time_type &t , time_type &dt )
	{
		return try_step_v1( system , x , t, dt );
	}

	/*
	 * Version 2 : try_step( sys , x , dxdt , t , dt )
	 *
	 * this version does not solve the forwarding problem, boost.range can not be used
	 */
	template< class System , class StateInOut , class DerivIn >
	controlled_step_result try_step( System system , StateInOut &x , const DerivIn &dxdt , time_type &t , time_type &dt )
	{
		m_xnew_resizer.adjust_size( x , boost::bind( &controlled_error_bs_type::resize_m_xnew< StateInOut > , boost::ref( *this ) , _1 ) );
		controlled_step_result res = try_step( system , x , dxdt , t , m_xnew.m_v , dt );
		if( ( res == success_step_size_increased ) || ( res == success_step_size_unchanged ) )
		{
			boost::numeric::odeint::copy( m_xnew.m_v , x );
		}
		return res;
	}

	/*
	 * Version 3 : try_step( sys , in , t , out , dt )
	 *
	 * this version does not solve the forwarding problem, boost.range can not be used
	 */
	template< class System , class StateIn , class StateOut >
	controlled_step_result try_step( System system , const StateIn &in , time_type &t , StateOut &out , time_type &dt )
	{
		typename boost::unwrap_reference< System >::type &sys = system;
		m_dxdt_resizer.adjust_size( in , boost::bind( &controlled_error_bs_type::resize_m_dxdt< StateIn > , boost::ref( *this ) , _1 ) );
		sys( in , m_dxdt.m_v , t );
		return try_step( system , in , m_dxdt.m_v , t , out , dt );
	}


    template< class System , class StateIn , class DerivIn , class StateOut >
	controlled_step_result try_step( System system , const StateIn &in , const DerivIn &dxdt , time_type &t , StateOut &out , time_type &dt )
	{
        typename boost::unwrap_reference< System >::type &sys = system;
		m_x_resizer.adjust_size( in , boost::bind( &controlled_error_bs_type::resize_m_x< StateIn > , boost::ref( *this ) , _1 ) );

        // actual algorithm ... 
        
        return step_size_decreased;
    }


    /* Resizer methods */

    template< class StateIn >
    bool resize_m_dxdt( const StateIn &x )
    {
        return adjust_size_by_resizeability( m_dxdt , x , typename wrapped_deriv_type::is_resizeable() );
    }

    template< class StateIn >
	bool resize_m_xnew( const StateIn &x )
	{
	    return adjust_size_by_resizeability( m_xnew , x , typename wrapped_state_type::is_resizeable() );
	}

    template< class StateIn >
	bool resize_m_x( const StateIn &x )
	{
        bool resized( false );
	    resized |= adjust_size_by_resizeability( m_xerr , x , typename wrapped_state_type::is_resizeable() );
        resized |= adjust_size_by_resizeability( m_xmp , x , typename wrapped_state_type::is_resizeable() );
        return resized;
	}

    template< class StateIn >
    bool resize_extrapolate( const StateIn &x )
    {
        bool resized( false );
        resized |= adjust_size_by_resizeability( m_c , x , typename wrapped_state_type::is_resizeable() );
    }

    template< class StateIn >
    void adjust_size( const StateIn &x )
    {
        resize_m_dxdt( x );
        resize_m_xnew( x );
        resize_m_x( x );
        m_midpoint.adjust_size();
    }


private:

    template< class System , class StateInOut >
	controlled_step_result try_step_v1( System system , StateInOut &x , time_type &t , time_type &dt )
	{
		typename boost::unwrap_reference< System >::type &sys = system;
		m_dxdt_resizer.adjust_size( x , boost::bind( &controlled_error_bs_type::resize_m_dxdt< StateInOut > , boost::ref( *this ) , _1 ) );
		sys( x , m_dxdt.m_v ,t );
		return try_step( system , x , m_dxdt.m_v , t , dt );
	}

    default_error_checker< value_type, algebra_type , operations_type > m_error_checker;
    explicit_midpoint< state_type , value_type , deriv_type , time_type , algebra_type , operations_type , resizer_type > m_midpoint;

    const unsigned short m_k_max;

    const time_type m_safety1;
    const time_type m_safety2;
    const time_type m_max_dt_factor;
    const time_type m_min_step_scale;
    const time_type m_max_step_scale;

    bool m_continuous_calls;
    bool m_decreased_step_during_last_call;

    time_type m_dt_last;
    time_type m_t_last;
    time_type m_current_eps;

    unsigned short m_current_k_max;
    unsigned short m_current_k_opt;

    resizer_type m_dxdt_resizer;
    resizer_type m_xnew_resizer;
    resizer_type m_x_resizer;

    wrapped_state_type m_xnew;
    wrapped_state_type m_xerr;
    wrapped_state_type m_xmp;
    wrapped_deriv_type m_dxdt;

    typedef std::vector< time_type > value_vector;
    typedef std::vector< value_vector > value_matrix;
    typedef std::vector< unsigned short > us_vector;

    value_vector m_error; // errors of repeated midpoint steps and extrapolations
    value_vector m_a; // stores the work (number of f calls) required for the orders
    value_matrix m_alpha; // stores convergence factor for stepsize adjustment
    us_vector m_interval_sequence;

    value_vector m_times;
    std::vector< wrapped_state_type > m_d;
    wrapped_state_type m_c;

};

}
}
}

#endif /* BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_ERROR_BS_HPP_ */
