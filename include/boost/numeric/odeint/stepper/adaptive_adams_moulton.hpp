#ifndef STEPPER_ADAMS_MOULTON_HPP_INCLUDED
#define STEPPER_ADAMS_MOULTON_HPP_INCLUDED

// helper class for controlled_adams_bashforth_moulton

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>
#include <boost/numeric/odeint/util/bind.hpp>

#include <boost/numeric/odeint/stepper/base/algebra_stepper_base.hpp>
#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/operations_dispatcher.hpp>

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
class adaptive_adams_moulton: public algebra_stepper_base< Algebra , Operations >
{
	public:
		static const size_t steps = Steps;
		typedef unsigned short order_type;
		static const order_type order = steps + 1;

		typedef State state_type;
		typedef Value value_type;
		typedef Deriv deriv_type;
		typedef Time time_type;
		typedef Resizer resizer_type;

		typedef Algebra algebra_type;
		typedef Operations operations_type;
		typedef algebra_stepper_base< Algebra , Operations > algebra_stepper_base_type;

		typedef state_wrapper<state_type> wrapped_state_type;

		typedef detail::adaptive_adams_coefficients<steps, deriv_type, time_type, algebra_type, operations_type> coeff_type;

		typedef adaptive_adams_moulton< Steps , State , Value , Deriv , Time , Resizer > stepper_type;

		adaptive_adams_moulton( const algebra_type &algebra = algebra_type())
		:algebra_stepper_base_type( algebra ), m_coeff()
		{};

		void do_step(coeff_type & coeff, const state_type & in,
					time_type t, state_type & out, time_type &dt)
		{
			coeff.poly.add_root(0);
			out = in;

			// integrating
			for(size_t i=0; i<coeff.m_effective_order; ++i)
			{
				time_type c = ((i!=coeff.m_effective_order-1)?coeff.m_c[i]:coeff.poly.evaluate_integrated(dt));
				this->m_algebra.for_each3(out, out, coeff.m_tss[i][coeff.m_effective_order-i-1].m_v, typename Operations::template scale_sum2<double, double>(1.0, c));
			}
			// std::cout << std::endl;
		};

	private:
		coeff_type m_coeff;
};

}
}
}

#endif