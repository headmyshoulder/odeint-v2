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

#include <boost/ref.hpp>

#include <boost/numeric/odeint/util/size_adjuster.hpp>
#include <boost/numeric/odeint/util/construct.hpp>
#include <boost/numeric/odeint/util/destruct.hpp>
#include <boost/numeric/odeint/util/copy.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>



namespace boost {
namespace numeric {
namespace odeint {


/*
 * Error checker for controlled_error_stepper
 *
 * ToDo: implement constructor with epsilons
 */
template
<
	class Value ,
	class Algebra = range_algebra ,
	class Operations = default_operations
>
class error_checker_standard
{
public:

	typedef Value value_type;
	typedef Algebra algebra_type;
	typedef Operations operations_type;


	error_checker_standard( void ) : m_eps_abs( 1E-6 ) , m_eps_rel( 1E-6 ) , m_a_x( 1.0 ) , m_a_dxdt( 1.0 )
	{}


	template< class State , class Deriv , class Err , class Time >
	value_type error( const State &x_old , const Deriv &dxdt_old , Err &x_err , const Time &dt )
	{
		// this overwrites x_err !
		algebra_type::for_each3( x_old , dxdt_old , x_err ,
					             typename operations_type::template rel_error< value_type >( m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt * detail::get_value( dt ) ) );

		value_type res = algebra_type::reduce( x_err , typename operations_type::template maximum< value_type >() , 0.0 );
		return res;
	}

private:

	value_type m_eps_abs;
	value_type m_eps_rel;
	value_type m_a_x;
	value_type m_a_dxdt;
};








/*
 * error stepper category dispatcher
 */
template<
    class ErrorStepper ,
    class ErrorChecker = error_checker_standard< typename ErrorStepper::value_type ,
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
class controlled_error_stepper< ErrorStepper , ErrorChecker , AdjustSizePolicy , explicit_error_stepper_tag >
{
	void initialize( void )
	{
		boost::numeric::odeint::construct( m_dxdt );
		boost::numeric::odeint::construct( m_xerr );
		boost::numeric::odeint::construct( m_xnew );
		m_dxdt_size_adjuster.register_state( 0 , m_dxdt );
		m_xerr_size_adjuster.register_state( 0 , m_xerr );
		m_xnew_size_adjuster.register_state( 0 , m_xnew );
	}

	void copy( const controlled_error_stepper &stepper )
	{
		boost::numeric::odeint::copy( stepper.m_dxdt , m_dxdt );
		boost::numeric::odeint::copy( stepper.m_xerr , m_xerr );
		boost::numeric::odeint::copy( stepper.m_xnew , m_xnew );
		m_max_rel_error = stepper.m_max_rel_error;
	}

public:

	typedef ErrorStepper stepper_type;
	typedef typename stepper_type::state_type state_type;
	typedef typename stepper_type::value_type value_type;
	typedef typename stepper_type::deriv_type deriv_type;
	typedef typename stepper_type::time_type time_type;
	typedef typename stepper_type::order_type order_type;
	typedef AdjustSizePolicy adjust_size_policy;
	typedef ErrorChecker error_checker_type;
	typedef controlled_stepper_tag stepper_category;



	controlled_error_stepper(
			const stepper_type &stepper = stepper_type() ,
			const error_checker_type &error_checker = error_checker_type()
			)
	: m_stepper( stepper ) , m_error_checker( error_checker ) ,
	  m_dxdt_size_adjuster() , m_xerr_size_adjuster() , m_xnew_size_adjuster() ,
	  m_dxdt() , m_xerr() , m_xnew() , m_max_rel_error()
	{
		initialize();
	}

	~controlled_error_stepper( void )
	{
		boost::numeric::odeint::destruct( m_dxdt );
		boost::numeric::odeint::destruct( m_xerr );
		boost::numeric::odeint::destruct( m_xnew );
	}

	controlled_error_stepper( const controlled_error_stepper &stepper )
	: m_stepper( stepper.m_stepper ) , m_error_checker( stepper.m_error_checker ) ,
	  m_dxdt_size_adjuster() , m_xerr_size_adjuster() , m_xnew_size_adjuster() ,
	  m_dxdt() , m_xerr() , m_xnew() , m_max_rel_error()
	{
		initialize();
		copy( stepper );
	}

	controlled_error_stepper& operator=( const controlled_error_stepper &stepper )
	{
		m_stepper = stepper.m_stepper;
		m_error_checker = stepper.m_error_checker;
		copy( stepper );
		return *this;
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
		m_xnew_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
		controlled_step_result res = try_step( system , x , dxdt , t , m_xnew , dt );
		if( ( res == success_step_size_increased ) || ( res == success_step_size_unchanged ) )
		{
			boost::numeric::odeint::copy( m_xnew , x );
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
		m_dxdt_size_adjuster.adjust_size_by_policy( in , adjust_size_policy() );
		sys( in , m_dxdt , t );
		return try_step( system , in , m_dxdt , t , out , dt );
	}


	/*
	 * Version 4 : try_step( sys , in , dxdt , t , out , dt )
	 *
	 * this version does not solve the forwarding problem, boost.range can not be used
	 */
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


	template< class System , class StateInOut >
	controlled_step_result try_step_v1( System system , StateInOut &x , time_type &t , time_type &dt )
	{
		typename boost::unwrap_reference< System >::type &sys = system;
		m_dxdt_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
		sys( x , m_dxdt ,t );
		return try_step( system , x , m_dxdt , t , dt );
	}


	stepper_type m_stepper;
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
class controlled_error_stepper< ErrorStepper , ErrorChecker , AdjustSizePolicy , explicit_error_stepper_fsal_tag >
{
	void initialize( void )
	{
		boost::numeric::odeint::construct( m_dxdt );
		boost::numeric::odeint::construct( m_xerr );
		boost::numeric::odeint::construct( m_xnew );
		boost::numeric::odeint::construct( m_dxdtnew );
		m_dxdt_size_adjuster.register_state( 0 , m_dxdt );
		m_xerr_size_adjuster.register_state( 0 , m_xerr );
		m_x_new_size_adjuster.register_state( 0 , m_xnew );
		m_dxdt_new_size_adjuster.register_state( 0 , m_dxdtnew );
	}

	void copy( const controlled_error_stepper &stepper )
	{
		boost::numeric::odeint::copy( stepper.m_dxdt , m_dxdt );
		boost::numeric::odeint::copy( stepper.m_xerr , m_xerr );
		boost::numeric::odeint::copy( stepper.m_xnew , m_xnew );
		boost::numeric::odeint::copy( stepper.m_dxdtnew , m_dxdtnew );
		m_first_call = stepper.m_first_call;
	}

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
            const stepper_type &stepper = stepper_type() ,
            const error_checker_type &error_checker = error_checker_type()
            )
    : m_stepper( stepper ) , m_error_checker( error_checker ) ,
      m_dxdt_size_adjuster() , m_xerr_size_adjuster() , m_x_new_size_adjuster() , m_dxdt_new_size_adjuster() ,
      m_dxdt() , m_xerr() , m_xnew() , m_dxdtnew() ,
      m_first_call( true )
    {
    	initialize();
    }

    controlled_error_stepper( const controlled_error_stepper &stepper )
    : m_stepper( stepper.m_stepper ) , m_error_checker( stepper.m_error_checker ) ,
      m_dxdt_size_adjuster() , m_xerr_size_adjuster() , m_x_new_size_adjuster() , m_dxdt_new_size_adjuster() ,
      m_dxdt() , m_xerr() , m_xnew() , m_dxdtnew() ,
      m_first_call( true )
    {
    	initialize();
    	copy( stepper );
    }

    ~controlled_error_stepper( void )
    {
        boost::numeric::odeint::destruct( m_dxdt );
        boost::numeric::odeint::destruct( m_xerr );
        boost::numeric::odeint::destruct( m_xnew );
        boost::numeric::odeint::destruct( m_dxdtnew );
    }

    controlled_error_stepper& operator=( const controlled_error_stepper &stepper )
    {
    	copy( stepper );
    	return *this;
    }





	/*
	 * Version 1 : try_step( sys , x , t , dt )
	 *
	 * The two overloads are needed in order to solve the forwarding problem
	 */
    template< class System , class StateInOut >
    controlled_step_result try_step( System system , StateInOut &x , time_type &t , time_type &dt )
    {
    	return try_step_v1( system , x , t , dt );
    }

    template< class System , class StateInOut >
    controlled_step_result try_step( System system , const StateInOut &x , time_type &t , time_type &dt )
    {
    	return try_step_v1( system , x , t , dt );
    }



	/*
	 * Version 2 : try_step( sys , in , t , out , dt );
	 *
	 * This version does not solve the forwarding problem, boost::range can not be used.
	 */
    template< class System , class StateIn , class StateOut >
    controlled_step_result try_step( System system , const StateIn &in , time_type &t , StateOut &out , time_type &dt )
    {
        if( m_dxdt_size_adjuster.adjust_size_by_policy( in , adjust_size_policy() ) || m_first_call )
        {
    		typename boost::unwrap_reference< System >::type &sys = system;
            sys( in , m_dxdt ,t );
            m_first_call = false;
        }
        return try_step( system , in , m_dxdt , t , out , dt );
    }


	/*
	 * Version 3 : try_step( sys , x , dxdt , t , dt )
	 *
	 * This version does not solve the forwarding problem, boost::range can not be used.
	 */
    template< class System , class StateInOut , class DerivInOut >
    controlled_step_result try_step( System system , StateInOut &x , DerivInOut &dxdt , time_type &t , time_type &dt )
    {
    	m_x_new_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
    	m_dxdt_new_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() );
    	controlled_step_result res = try_step( system , x , dxdt , t , m_xnew , m_dxdtnew , dt );
    	if( ( res == success_step_size_increased ) || ( res == success_step_size_unchanged) )
    	{
    		boost::numeric::odeint::copy( m_xnew , x );
    		boost::numeric::odeint::copy( m_dxdtnew , dxdt );
    	}
    	return res;
    }


	/*
	 * Version 3 : try_step( sys , in , dxdt , t , out , dt )
	 *
	 * This version does not solve the forwarding problem, boost::range can not be used.
	 */
    template< class System , class StateIn , class DerivIn , class StateOut , class DerivOut >
    controlled_step_result try_step( System system , const StateIn &in , const DerivIn &dxdt_in , time_type &t ,
    		StateOut &out , DerivOut &dxdt_out , time_type &dt )
    {
        using std::max;
        using std::min;
        using std::pow;

        m_xerr_size_adjuster.adjust_size_by_policy( in , adjust_size_policy() );

        //fsal: m_stepper.get_dxdt( dxdt );
        //fsal: m_stepper.do_step( sys , x , dxdt , t , dt , m_x_err );
        m_stepper.do_step( system , in , dxdt_in , t , out , dxdt_out , dt , m_xerr );

        // this potentially overwrites m_x_err! (standard_error_checker does, at least)
        value_type max_rel_err = m_error_checker.error( in , dxdt_in , m_xerr , dt );

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
        changed |= m_x_new_size_adjuster.adjust_size( x );
        changed |= m_dxdt_new_size_adjuster.adjust_size( x );
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

    template< class System , class StateInOut >
    controlled_step_result try_step_v1( System system , StateInOut &x , time_type &t , time_type &dt )
    {
        if( m_dxdt_size_adjuster.adjust_size_by_policy( x , adjust_size_policy() ) || m_first_call )
        {
    		typename boost::unwrap_reference< System >::type &sys = system;
            sys( x , m_dxdt ,t );
            m_first_call = false;
        }
        return try_step( system , x , m_dxdt , t , dt );
    }


    stepper_type m_stepper;
    error_checker_type m_error_checker;

    size_adjuster< deriv_type , 1 > m_dxdt_size_adjuster;
    size_adjuster< state_type , 1 > m_xerr_size_adjuster;
    size_adjuster< state_type , 1 > m_x_new_size_adjuster;
    size_adjuster< deriv_type , 1 > m_dxdt_new_size_adjuster;

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
