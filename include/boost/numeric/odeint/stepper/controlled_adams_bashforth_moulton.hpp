#ifndef STEPPER_CONTROLLED_ADAMS_BASHFORTH_MOULTON_HPP_INCLUDED
#define STEPPER_CONTROLLED_ADAMS_BASHFORTH_MOULTON_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>

#include <boost/numeric/odeint/stepper/adaptive_adams_bashforth_moulton.hpp>
#include <boost/numeric/odeint/stepper/detail/pid_step_adjuster.hpp>

#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>
#include <boost/numeric/odeint/util/bind.hpp>

namespace boost{
namespace numeric{
namespace odeint {

// template<
// class ErrorStepper,
// class ErrorChecker,
// class StepAdjuster,
// class Resizer
// >
template<
class ErrorStepper,
class StepAdjuster = detail::pid_step_adjuster<ErrorStepper::order_value, typename ErrorStepper::state_type, typename ErrorStepper::time_type>,
class Resizer = initially_resizer
>
class controlled_adams_bashforth_moulton
{
	public:
		typedef ErrorStepper stepper_type;
		typedef typename stepper_type::state_type state_type;
		typedef typename stepper_type::value_type value_type;
		typedef typename stepper_type::deriv_type deriv_type;
		typedef typename stepper_type::time_type time_type;

		typedef typename stepper_type::algebra_type algebra_type;
		typedef typename stepper_type::operations_type operations_type;
		typedef Resizer resizer_type;

		// typedef ErrorChecker error_checker_type;
		typedef StepAdjuster step_adjuster_type;
		typedef controlled_stepper_tag stepper_category;

		typedef typename stepper_type::wrapped_state_type wrapped_state_type;
		typedef typename stepper_type::wrapped_deriv_type wrapped_deriv_type;

		typedef controlled_adams_bashforth_moulton<ErrorStepper, StepAdjuster, Resizer> controlled_stepper_type;

		controlled_adams_bashforth_moulton(
			step_adjuster_type step_adjuster = step_adjuster_type()
			)
		:m_stepper(), m_coeff(m_stepper.coeff()), 
		m_dxdt_resizer(), m_xerr_resizer(), m_xnew_resizer(),
		m_step_adjuster(step_adjuster)
		{};

		template<class System>
		controlled_step_result try_step(System system, state_type & inOut, time_type &t, time_type &dt)
		{
			m_xnew_resizer.adjust_size( inOut , detail::bind( &controlled_stepper_type::template resize_xnew_impl< state_type > , detail::ref( *this ) , detail::_1 ) );

			controlled_step_result res = try_step(system, inOut, t, m_xnew.m_v, dt);

			if(res == success)
				boost::numeric::odeint::copy( m_xnew.m_v , inOut);

			return res;
		};

		template<class System>
		controlled_step_result try_step(System system, const state_type & in, time_type &t, state_type & out, time_type &dt)
		{
			m_xerr_resizer.adjust_size( in , detail::bind( &controlled_stepper_type::template resize_xerr_impl< state_type > , detail::ref( *this ) , detail::_1 ) );
			m_stepper.do_step_impl(system, m_coeff, in, t, out, dt, m_xerr.m_v);

			double ratio = m_step_adjuster.adjust_stepsize(m_xerr.m_v, dt);

			if(ratio >= 0.9)
			{
				m_dxdt_resizer.adjust_size( in , detail::bind( &controlled_stepper_type::template resize_dxdt_impl< state_type > , detail::ref( *this ) , detail::_1 ) );

				system(out, m_dxdt.m_v, t+dt);
				m_coeff.step(m_dxdt.m_v, t+dt);
				m_coeff.confirm();

				t += dt;
				dt *= ratio;

				return success;
			}
			else
			{
				dt *= ratio;
				return fail;
			}
		};

	private:
		template< class StateType >
		bool resize_dxdt_impl( const StateType &x )
		{
			return adjust_size_by_resizeability( m_dxdt, x, typename is_resizeable<deriv_type>::type() );
		};
		template< class StateType >
		bool resize_xerr_impl( const StateType &x )
		{
			return adjust_size_by_resizeability( m_xerr, x, typename is_resizeable<state_type>::type() );
		};
		template< class StateType >
		bool resize_xnew_impl( const StateType &x )
		{
			return adjust_size_by_resizeability( m_xnew, x, typename is_resizeable<state_type>::type() );
		};

		stepper_type m_stepper;
		typename stepper_type::coeff_type &m_coeff;

		wrapped_deriv_type m_dxdt;
		wrapped_state_type m_xerr;
		wrapped_state_type m_xnew;

		resizer_type m_dxdt_resizer;
		resizer_type m_xerr_resizer;
		resizer_type m_xnew_resizer;

		// error_checker_type m_error_checker;
		step_adjuster_type m_step_adjuster;
};

}
}
}

#endif