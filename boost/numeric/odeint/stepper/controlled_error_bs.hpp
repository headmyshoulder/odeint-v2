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
            time_type abs_err, time_type rel_err,
            time_type factor_x, time_type factor_dxdt )
        : m_error_checker( abs_err, rel_err, factor_x, factor_dxdt ),
          m_k_max(8),
          m_safety1(0.25), m_safety2(0.7),
          m_max_dt_factor( 0.1 ), m_min_step_scale(5E-5), m_max_step_scale(0.7),
          m_continuous_calls(false), m_decreased_step_during_last_call( false ),
          m_dt_last( 1.0E30),
          m_current_eps( -1.0 )
    {
        m_error.resize(m_k_max);
        m_a.resize(m_k_max+1);
        m_alpha.resize(m_k_max); // k_max * k_max matrix
        typename value_matrix::iterator it = m_alpha.begin();
        while( it != m_alpha.end() )
            (*it++).resize(m_k_max);
        m_interval_sequence.resize(m_k_max+1);
        for( unsigned short i = 1; i <= m_k_max+1; i++ )
            m_interval_sequence[i-1] = (2*i);

        m_times.resize(m_k_max);
        m_d.resize(m_k_max);
    }

    template< class System , class StateIn , class DerivIn , class StateOut >
    controlled_step_result try_step( System system , const StateIn &in , const DerivIn &dxdt , time_type &t , StateOut &out , time_type &dt )
    {

    }

private:

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

    wrapped_state_type m_x0;
    wrapped_state_type m_xerr;
    wrapped_state_type m_x_mp;
    wrapped_state_type m_x_scale;
    wrapped_deriv_type m_dxdt;

    typedef std::vector< time_type > value_vector;
    typedef std::vector< std::vector< time_type > > value_matrix;
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
