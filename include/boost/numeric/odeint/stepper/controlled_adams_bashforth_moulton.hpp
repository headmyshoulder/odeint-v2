#ifndef STEPPER_CONTROLLED_ADAMS_BASHFORTH_MOULTON_HPP_INCLUDED
#define STEPPER_CONTROLLED_ADAMS_BASHFORTH_MOULTON_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/adaptive_adams_bashforth.hpp>
#include <boost/numeric/odeint/stepper/adaptive_adams_moulton.hpp>
#include <boost/numeric/odeint/stepper/detail/adaptive_adams_coefficients.hpp>
#include <boost/numeric/odeint/stepper/detail/pid_step_adjuster.hpp>

#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/operations_dispatcher.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/util/bind.hpp>

#include <math.h>

namespace boost {
namespace numeric {
namespace odeint {



template<
size_t Steps,
class State,
class Value,
class Deriv = State,
class Time = Value,
class StepAdjuster = detail::pid_step_adjuster<Steps, State, Time>,
class Algebra = typename algebra_dispatcher< State >::algebra_type,
class Operations = typename operations_dispatcher< State >::operations_type ,
class Resizer = initially_resizer
>
class controlled_adams_bashforth_moulton
{
	public:
		static const size_t steps = Steps;
		typedef unsigned short order_type;
		static const order_type order = steps;

		typedef State state_type;
		typedef Value value_type;
		typedef Deriv deriv_type;
		typedef Time time_type;
		typedef StepAdjuster step_adjuster_type;
		typedef Resizer resizer_type;

		typedef Algebra algebra_type;
		typedef Operations operations_type;

		typedef state_wrapper<state_type> wrapped_state_type;
		typedef state_wrapper<deriv_type> wrapped_deriv_type;

		typedef detail::adaptive_adams_coefficients<order, deriv_type, time_type> coeff_type;

		typedef detail::rotating_buffer<state_type, steps> error_storage_type;

		typedef adaptive_adams_bashforth<order, state_type, value_type> aab_type;
		typedef adaptive_adams_moulton<order, state_type, value_type> aam_type;
		typedef controlled_adams_bashforth_moulton< Steps , State , Value , Deriv , Time, StepAdjuster, Algebra, Operations, Resizer > stepper_type;

		typedef explicit_controlled_stepper_tag stepper_category;

		controlled_adams_bashforth_moulton()
		:m_aab(), m_aam(m_aab.algebra()), m_coeff(),
		m_xerr_resizer(), m_dxdt_resizer(), m_xnew_resizer(),
		m_step_adjuster()
		{};

		controlled_adams_bashforth_moulton(const algebra_type &algebra)
		:m_aab(algebra), m_aam(m_aab.algebra()), m_coeff(),
		m_xerr_resizer(), m_dxdt_resizer(), m_xnew_resizer(),
		m_step_adjuster()
		{};

		template<class System>
		controlled_step_result try_step(System system, state_type & inOut, time_type &t, time_type &dt)
		{
			m_xnew_resizer.adjust_size( inOut , detail::bind( &stepper_type::template resize_xnew_impl< state_type > , detail::ref( *this ) , detail::_1 ) );

			controlled_step_result res = try_step(system, inOut, t, m_xnew.m_v, dt);

			if(res == success)
				boost::numeric::odeint::copy( m_xnew.m_v , inOut);

			return res;
		};

		template<class System>
		controlled_step_result try_step(System system, const state_type & in, time_type &t, state_type & out, time_type &dt)
		{
			m_dxdt_resizer.adjust_size( in , detail::bind( &stepper_type::template resize_dxdt_impl< state_type > , detail::ref( *this ) , detail::_1 ) );

			if(m_coeff.m_effective_order == 1)
			{
				system(in, m_dxdt.m_v, t);

				m_coeff.step(m_dxdt.m_v, t);
				m_coeff.confirm();
			}
			// predict
			m_aab.do_step_impl(m_coeff, in, t, out, dt);

			// evaluate
			system(out, m_dxdt.m_v, t + dt);
			m_coeff.step(m_dxdt.m_v, t + dt);

			m_xerr_resizer.adjust_size( in , detail::bind( &stepper_type::template resize_xerr_impl< state_type > , detail::ref( *this ) , detail::_1 ) );
			boost::numeric::odeint::copy( out, m_xerr.m_v);

			// correct
			m_aam.do_step(m_coeff, in, t, out, dt);
			m_aab.algebra().for_each3(m_xerr.m_v, m_xerr.m_v, out, typename Operations::template scale_sum2<double, double>(1.0, -1.0));

			double ratio = m_step_adjuster.adjust_stepsize(m_xerr.m_v, dt);

			if(ratio >= 0.9)
			{
				// evaluate
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
			return adjust_size_by_resizeability( m_xerr, x, typename is_resizeable<deriv_type>::type() );
		};
		template< class StateType >
		bool resize_xnew_impl( const StateType &x )
		{
			return adjust_size_by_resizeability( m_xnew, x, typename is_resizeable<state_type>::type() );
		};

		aab_type m_aab;
		aam_type m_aam;

		resizer_type m_dxdt_resizer;
		resizer_type m_xerr_resizer;
		resizer_type m_xnew_resizer;

		wrapped_deriv_type m_dxdt;
		wrapped_state_type m_xnew;
		wrapped_state_type m_xerr;

		coeff_type m_coeff;

		step_adjuster_type m_step_adjuster;
};

}}}

#endif