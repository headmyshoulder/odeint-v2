/*
 boost header: NUMERIC_ODEINT_STEPPER/controlled_error_stepper.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_ERROR_STEPPER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_ERROR_STEPPER_HPP_INCLUDED

#include <cmath>

#include <boost/noncopyable.hpp>
#include <boost/ref.hpp>

#include <boost/numeric/odeint/stepper/adjust_size.hpp>
#include <boost/numeric/odeint/stepper/error_checker.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

namespace boost {
namespace numeric {
namespace odeint {

typedef enum
{
    success_step_size_unchanged ,
    step_size_decreased ,
    success_step_size_increased
} controlled_step_result;


/*
 * error stepper category dispatcher
 */
template<
    class ErrorStepper ,
    class ErrorChecker = error_checker_standard< typename ErrorStepper::state_type ,
                                                   typename ErrorStepper::time_type ,
                                                   typename ErrorStepper::algebra_type ,
                                                   typename ErrorStepper::operations_type > ,
    class AdjustSizePolicy = typename ErrorStepper::adjust_size_policy ,
    class ErrorStepperCategory = typename ErrorStepper::stepper_category
>
class controlled_error_stepper { };




/*
 * explicit stepper version
 */
template<
	class ErrorStepper ,
	class ErrorChecker ,
	class AdjustSizePolicy
	>
class controlled_error_stepper< ErrorStepper , ErrorChecker , AdjustSizePolicy , explicit_error_stepper_tag > : boost::noncopyable
{
public:

	typedef ErrorStepper stepper_type;
	typedef typename stepper_type::state_type state_type;
	typedef typename stepper_type::value_type value_type;
	typedef typename stepper_type::deriv_type deriv_type;
	typedef typename stepper_type::time_type time_type;
	typedef typename stepper_type::order_type order_type;
	typedef AdjustSizePolicy adjust_size_policy;
	typedef ErrorChecker error_checker_type;



	controlled_error_stepper(
			stepper_type &stepper ,
			const error_checker_type &error_checker = error_checker_type()
			)
	: m_stepper( stepper ) , m_error_checker( error_checker ) ,
	  m_dxdt_size_adjuster() , m_xerr_size_adjuster() , m_xnew_size_adjuster() ,
	  m_dxdt() , m_xerr() , m_xnew() , m_max_rel_error()
	{
		boost::numeric::odeint::construct( m_dxdt );
		boost::numeric::odeint::construct( m_xerr );
		boost::numeric::odeint::construct( m_xnew );
		m_dxdt_size_adjuster.register_state( 0 , m_dxdt );
		m_xerr_size_adjuster.register_state( 0 , m_xerr );
		m_xnew_size_adjuster.register_state( 0 , m_xnew );
	}

	~controlled_error_stepper( void )
	{
		boost::numeric::odeint::destruct( m_dxdt );
		boost::numeric::odeint::destruct( m_xerr );
		boost::numeric::odeint::destruct( m_xnew );
	}




	// try_step( sys , x , t , dt )
	template< class System , class StateInOut >
	controlled_step_result try_step( System system , StateInOut &x , time_type &t , time_type &dt )
	{
		typename boost::unwrap_reference< System >::type &sys = system;
		m_dxdt_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
		sys( x , m_dxdt ,t );
		return try_step( system , x , m_dxdt , t , dt );
	}

	// try_step( sys , x , dxdt , t , dt )
	template< class System , class StateInOut , class DerivIn >
	controlled_step_result try_step( System system , StateInOut &x , const DerivIn &dxdt , time_type &t , time_type &dt )
	{
		m_xnew_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
		controlled_step_result res = try_step( system , x , dxdt , t , m_xnew , dt );
		if( ( res == success_step_size_increased ) || ( res == success_step_size_unchanged ) )
		{
			boost::numeric::odeint::copy( m_xnew , x );
		}
		return res;
	}

	// try_step( sys , in , t , out , dt )
	template< class System , class StateIn , class StateOut >
	controlled_step_result try_step( System system , const StateIn &in , time_type &t , StateOut &out , time_type &dt )
	{
		typename boost::unwrap_reference< System >::type &sys = system;
		m_dxdt_size_adjuster.adjust_size_by_policy( in , adjust_size_policy() );
		system( in , m_dxdt , t );
		return try_step( system , in , m_dxdt , t , out , dt );
	}


	// try_step( sys , in , dxdt , t , out , dt )
	template< class System , class StateIn , class DerivIn , class StateOut >
	controlled_step_result try_step( System system , const StateIn &in , const DerivIn &dxdt , time_type &t , StateOut &out , time_type &dt )
	{
		using std::max;
		using std::min;
		using std::pow;

		m_xerr_size_adjuster.adjust_size_by_policy( in , adjust_size_policy() );

		// do one step with error calculation
		m_stepper.do_step( system , in , dxdt , t , out , dt , m_xerr );

		m_max_rel_error = m_error_checker.error( in , dxdt , m_xerr , dt );

		if( m_max_rel_error > 1.1 )
		{
			// error too large - decrease dt ,limit scaling factor to 0.2 and reset state
			dt *= max( 0.9 * pow( m_max_rel_error , -1.0 / ( m_stepper.error_order() - 1.0 ) ) , 0.2 );
			return step_size_decreased;
		}
		else
		{
			if( m_max_rel_error < 0.5 )
			{
				//error too small - increase dt and keep the evolution and limit scaling factor to 5.0
				t += dt;
				dt *= min( 0.9 * pow( m_max_rel_error , -1.0 / m_stepper.stepper_order() ) , 5.0 );
				return success_step_size_increased;
			}
			else
			{
				t += dt;
				return success_step_size_unchanged;
			}
		}
	}

	value_type last_error( void ) const
	{
		return m_max_rel_error;
	}



	template< class StateType >
	void adjust_size( const StateType &x )
	{
        m_dxdt_size_adjuster.adjust_size( x );
        m_xerr_size_adjuster.adjust_size( x );
        m_xnew_size_adjuster.adjust_size( x );
        m_stepper.adjust_size( x );
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

	stepper_type &m_stepper;
	error_checker_type m_error_checker;

	size_adjuster< deriv_type , 1 > m_dxdt_size_adjuster;
	size_adjuster< state_type , 1 > m_xerr_size_adjuster;
	size_adjuster< state_type , 1 > m_xnew_size_adjuster;

	deriv_type m_dxdt;
	state_type m_xerr;
	state_type m_xnew;
	value_type m_max_rel_error;
};










/*
 * explicit stepper fsal version
 *
 * ToDo : introduce the same functions as for the above stepper
 */
template<
    class ErrorStepper ,
    class ErrorChecker ,
	class AdjustSizePolicy
    >
class controlled_error_stepper< ErrorStepper , ErrorChecker , AdjustSizePolicy , explicit_error_stepper_fsal_tag > : boost::noncopyable
{
public:

    typedef ErrorStepper stepper_type;
    typedef typename stepper_type::state_type state_type;
    typedef typename stepper_type::value_type value_type;
    typedef typename stepper_type::deriv_type deriv_type;
    typedef typename stepper_type::time_type time_type;
    typedef typename stepper_type::order_type order_type;
    typedef AdjustSizePolicy adjust_size_policy;
    typedef ErrorChecker error_checker_type;

    controlled_error_stepper(
            stepper_type &stepper ,
            const error_checker_type &error_checker = error_checker_type()
            )
    : m_stepper( stepper ) , m_error_checker( error_checker ) ,
      m_dxdt_size_adjuster() , m_xerr_size_adjuster() , m_new_size_adjuster() ,
      m_dxdt() , m_xerr() , m_xnew() , m_dxdtnew() ,
      m_first_call( true )
    {
        boost::numeric::odeint::construct( m_dxdt );
        boost::numeric::odeint::construct( m_xerr );
        boost::numeric::odeint::construct( m_xnew );
        boost::numeric::odeint::construct( m_dxdtnew );
        m_dxdt_size_adjuster.register_state( 0 , m_dxdt );
        m_xerr_size_adjuster.register_state( 0 , m_xerr );
        m_new_size_adjuster.register_state( 0 , m_xnew );
        m_new_size_adjuster.register_state( 1 , m_dxdtnew );
    }

    ~controlled_error_stepper( void )
    {
        boost::numeric::odeint::destruct( m_dxdt );
        boost::numeric::odeint::destruct( m_xerr );
        boost::numeric::odeint::destruct( m_xnew );
        boost::numeric::odeint::destruct( m_dxdtnew );
    }





    // try_step( sys , x , t , dt )
    template< class System >
    controlled_step_result try_step( System system , state_type &x , time_type &t , time_type &dt )
    {
        if( m_dxdt_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() ) || m_first_call )
        {
    		typename boost::unwrap_reference< System >::type &sys = system;
            sys( x , m_dxdt ,t );
            m_first_call = false;
        }
        return try_step( system , x , m_dxdt , t , dt );
    }

    // try_step( sys , in , t , out , dt );
    template< class System >
    controlled_step_result try_step( System system , const state_type &in , time_type &t , state_type &out , time_type &dt )
    {
        if( m_dxdt_size_adjuster.adjust_size_by_policy( in , adjust_size_policy() ) || m_first_call )
        {
    		typename boost::unwrap_reference< System >::type &sys = system;
            sys( in , m_dxdt ,t );
            m_first_call = false;
        }
        return try_step( system , in , m_dxdt , t , out , dt );
    }

    // try_step( sys , x , dxdt , t , dt )
    template< class System >
    controlled_step_result try_step( System system , state_type &x , state_type &dxdt , time_type &t , time_type &dt )
    {
    	m_new_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
    	controlled_step_result res = try_step( system , x , dxdt , t , m_xnew , m_dxdtnew , dt );
    	if( ( res == success_step_size_increased ) || ( res == success_step_size_unchanged) )
    	{
    		boost::numeric::odeint::copy( m_xnew , x );
    		boost::numeric::odeint::copy( m_dxdtnew , dxdt );
    	}
    	return res;
    }

    // try_step( sys , in , dxdt , t , out , dt )
    template< class System >
    controlled_step_result try_step( System system , const state_type &in , const state_type &dxdt_in , time_type &t ,
    		state_type &out , state_type &dxdt_out , time_type &dt )
    {
        using std::max;
        using std::min;
        using std::pow;

        m_xerr_size_adjuster.adjust_size_by_policy( in , adjust_size_policy() );

        //fsal: m_stepper.get_dxdt( dxdt );
        //fsal: m_stepper.do_step( sys , x , dxdt , t , dt , m_x_err );
        m_stepper.do_step( system , in , dxdt_in , t , out , dxdt_out , dt , m_xerr );

        // this potentially overwrites m_x_err! (standard_error_checker does, at least)
        time_type max_rel_err = m_error_checker.error( in , dxdt_in , m_xerr , dt );

        if( max_rel_err > 1.1 )
        {
            // error too large - decrease dt ,limit scaling factor to 0.2 and reset state
            dt *= max( 0.9 * pow( max_rel_err , -1.0 / ( m_stepper.error_order() - 1.0 ) ) , 0.2 );
            return step_size_decreased;
        }
        else
        {
            if( max_rel_err < 0.5 )
            {
                //error too small - increase dt and keep the evolution and limit scaling factor to 5.0
                t += dt;
                dt *= min( 0.9 * pow( max_rel_err , -1.0 / m_stepper.stepper_order() ) , 5.0 );
                return success_step_size_increased;
            }
            else
            {
                t += dt;
                return success_step_size_unchanged;
            }
        }
    }




    template< class StateType >
    void adjust_size( const StateType &x )
    {
        bool changed = false;
        changed |= m_dxdt_size_adjuster.adjust_size( x );
        changed |= m_xerr_size_adjuster.adjust_size( x );
        changed |= m_stepper.adjust_size( x );
        if( changed )
            m_first_call = true;
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

    stepper_type &m_stepper;
    error_checker_type m_error_checker;

    size_adjuster< deriv_type , 1 > m_dxdt_size_adjuster;
    size_adjuster< state_type , 1 > m_xerr_size_adjuster;
    size_adjuster< state_type , 2 > m_new_size_adjuster;

    deriv_type m_dxdt;
    state_type m_xerr;
    state_type m_xnew;
    deriv_type m_dxdtnew;
    bool m_first_call;
};




} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_ERROR_STEPPER_HPP_INCLUDED
