#ifndef STEPPER_ADAMS_BASHFORTH_HPP_INCLUDED
#define STEPPER_ADAMS_BASHFORTH_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/detail/polynomial.hpp>
#include <boost/numeric/odeint/stepper/detail/adaptive_adams_coefficients.hpp>

#include <boost/numeric/odeint/stepper/base/algebra_stepper_base.hpp>
#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/operations_dispatcher.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

#include <boost/numeric/odeint/util/copy.hpp>
#include <boost/numeric/odeint/util/bind.hpp>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

#include <iostream>

namespace boost {
namespace numeric {
namespace odeint {

template<
size_t Steps,
class State,
class Value = double,
class Deriv = State,
class Time = Value,
class Algebra = typename algebra_dispatcher< State >::algebra_type,
class Operations = typename operations_dispatcher< State >::operations_type ,
class Resizer = initially_resizer
>
class adaptive_adams_bashforth: public algebra_stepper_base< Algebra , Operations >
{
	public:
		static const size_t steps = Steps;
		typedef unsigned short order_type;
		static const order_type order_value = steps;

		typedef State state_type;
		typedef Value value_type;
		typedef Deriv deriv_type;
		typedef Time time_type;
		typedef Resizer resizer_type;

		typedef Algebra algebra_type;
		typedef Operations operations_type;
		typedef algebra_stepper_base< Algebra , Operations > algebra_stepper_base_type;

		typedef state_wrapper<state_type> wrapped_state_type;
		typedef state_wrapper<deriv_type> wrapped_deriv_type;
		typedef stepper_tag stepper_category;

		typedef detail::adaptive_adams_coefficients<Steps, deriv_type, time_type, algebra_type, operations_type> coeff_type;

		typedef adaptive_adams_bashforth< Steps , State , Value , Deriv , Time , Algebra, Operations, Resizer > stepper_type;

		adaptive_adams_bashforth( const algebra_type &algebra = algebra_type() )
		:algebra_stepper_base_type( algebra ) ,
		m_dxdt_resizer(), m_xnew_resizer()
		{};

		order_type order() const
		{
			return order_value;
		};

		template<class System>
		void do_step(System system, state_type & inOut, time_type t, time_type dt)
		{
			m_xnew_resizer.adjust_size( inOut , detail::bind( &stepper_type::template resize_xnew_impl< state_type > , detail::ref( *this ) , detail::_1 ) );

			do_step(system, inOut, t, m_xnew.m_v, dt);
			boost::numeric::odeint::copy( m_xnew.m_v , inOut);
		};

		template<class System>
		void do_step(System system, const state_type & in, time_type t, state_type & out, time_type dt)
		{
			m_dxdt_resizer.adjust_size( in , detail::bind( &stepper_type::template resize_dxdt_impl< state_type > , detail::ref( *this ) , detail::_1 ) );

			system(in, m_dxdt.m_v, t);
			m_coeff.step(m_dxdt.m_v, t);
			m_coeff.confirm();

			do_step_impl(m_coeff, in, t, out, dt);
		};

		void do_step_impl(coeff_type & coeff, const state_type & in, time_type t, state_type & out, time_type dt)
		{
			coeff.poly.reset();
			out = in;

			// integrating
			for(size_t i=0; i<coeff.m_effective_order-1; ++i)
			{
				if(i>0)
				{
					coeff.poly.add_root(coeff.m_ts[coeff.m_effective_order -1 -i] - coeff.m_ts[0]);
				}

				time_type c = coeff.poly.evaluate_integrated(dt);
				coeff.m_c[i] = c;

				// predict next state
				this->m_algebra.for_each3(out, out, coeff.m_ss[i][coeff.m_effective_order-i-2].m_v, typename Operations::template scale_sum2<double, double>(1.0, c));
			}
		};

		const coeff_type& coeff() const
		{
			return m_coeff;
		};

		coeff_type& coeff()
		{
			return m_coeff;
		};

		void reset()
		{
			m_coeff.reset();
		};

	private:
		template< class StateType >
		bool resize_dxdt_impl( const StateType &x )
		{
			return adjust_size_by_resizeability( m_dxdt, x, typename is_resizeable<deriv_type>::type() );
		};
		template< class StateType >
		bool resize_xnew_impl( const StateType &x )
		{
			return adjust_size_by_resizeability( m_xnew, x, typename is_resizeable<state_type>::type() );
		};

		resizer_type m_dxdt_resizer;
		resizer_type m_xnew_resizer;

		coeff_type m_coeff;
		wrapped_deriv_type m_dxdt;
		wrapped_state_type m_xnew;
};

}
}
}

#endif
