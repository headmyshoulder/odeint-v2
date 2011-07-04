/*
 * explicit_error_generic_rk.hpp
 *
 *  Created on: Jun 3, 2011
 *      Author: mario
 */

#ifndef EXPLICIT_ERROR_GENERIC_RK_HPP_
#define EXPLICIT_ERROR_GENERIC_RK_HPP_

#include <boost/numeric/odeint/stepper/base/explicit_stepper_and_error_stepper_base.hpp>

#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/stepper/detail/generic_rk_algorithm.hpp>
#include <boost/numeric/odeint/stepper/detail/generic_rk_call_algebra.hpp>
#include <boost/numeric/odeint/stepper/detail/generic_rk_operations.hpp>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;


namespace boost {
namespace numeric {
namespace odeint {

template<
    size_t StageCount,
    size_t Order,
    size_t StepperOrder ,
    size_t ErrorOrder ,
    class State ,
    class Value = double ,
    class Deriv = State ,
    class Time = Value ,
    class Algebra = range_algebra ,
    class Operations = default_operations ,
    class AdjustSizePolicy = adjust_size_initially_tag
    >
class explicit_error_generic_rk
: public explicit_stepper_and_error_stepper_base<
      explicit_error_generic_rk< StageCount , Order , StepperOrder , ErrorOrder , State ,
                                 Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > ,
      Order , StepperOrder , ErrorOrder , State , Value , Deriv , Time , Algebra ,
      Operations , AdjustSizePolicy >
{

public:

    typedef explicit_stepper_and_error_stepper_base<
            explicit_error_generic_rk< StageCount , Order , StepperOrder , ErrorOrder , State ,
                                       Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > ,
            Order , StepperOrder , ErrorOrder , State , Value , Deriv , Time , Algebra ,
            Operations , AdjustSizePolicy > stepper_base_type;

    typedef typename stepper_base_type::state_type state_type;
    typedef typename stepper_base_type::value_type value_type;
    typedef typename stepper_base_type::deriv_type deriv_type;
    typedef typename stepper_base_type::time_type time_type;
    typedef typename stepper_base_type::algebra_type algebra_type;
    typedef typename stepper_base_type::operations_type operations_type;
    typedef typename stepper_base_type::adjust_size_policy adjust_size_policy;
    typedef typename stepper_base_type::stepper_type stepper_type;

    typedef detail::generic_rk_algorithm< StageCount , Value , Algebra , Operations > rk_algorithm_type;

    typedef typename rk_algorithm_type::coef_a_type coef_a_type;
    typedef typename rk_algorithm_type::coef_b_type coef_b_type;
    typedef typename rk_algorithm_type::coef_c_type coef_c_type;

    static const size_t stage_count = StageCount;

private:

    void initialize( void )
    {
        boost::numeric::odeint::construct( m_x_tmp );
        m_state_adjuster.register_state( 0 , m_x_tmp );
        for( size_t i = 0 ; i < StageCount-1 ; ++i )
        {
            boost::numeric::odeint::construct( m_F[i] );
            m_deriv_adjuster.register_state( i , m_F[i] );
        }
    }

    void copy( const explicit_error_generic_rk &rk )
    {
        boost::numeric::odeint::copy( rk.m_x_tmp , m_x_tmp );
        for( size_t i = 0 ; i < StageCount-1 ; ++i )
        {
            boost::numeric::odeint::copy( rk.m_F[i] , m_F[i] );
        }
    }

public:

    // we use an explicit_generic_rk to do the normal rk step
    // and add a separate calculation of the error estimate afterwards
    explicit_error_generic_rk( const coef_a_type &a ,
                                  const coef_b_type &b ,
                                  const coef_b_type &b2 ,
                                  const coef_c_type &c )
        : m_rk_algorithm( a , b , c ) , m_b2( b2 ) , m_x_tmp()
    {
        initialize();
    }

    explicit_error_generic_rk( const explicit_error_generic_rk &rk )
        : stepper_base_type( rk )  , m_rk_algorithm( rk.m_rk_algorithm ) , m_b2( rk.m_b2 ) , m_x_tmp()
    {
        initialize();
        copy( rk );
    }

    explicit_error_generic_rk& operator=( const explicit_error_generic_rk &rk )
    {
        stepper_base_type::operator=( rk );
        copy( rk );
        return *this;
    }

    ~explicit_error_generic_rk( void )
   {
       boost::numeric::odeint::destruct( m_x_tmp );
       for( size_t i = 0 ; i < StageCount-1 ; ++i )
       {
           boost::numeric::odeint::destruct( m_F[i] );
       }
   }

    template< class System , class StateIn , class DerivIn , class StateOut , class Err >
    void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt ,
            const time_type &t , StateOut &out , const time_type &dt , Err &xerr )
    {
        // normal step
        do_step_impl( system , in , dxdt , t , out , dt );

        // additionally, perform the error calculation
        detail::template generic_rk_call_algebra< StageCount , algebra_type >()( xerr , dxdt , m_F ,
                    detail::generic_rk_scale_sum_err< StageCount , operations_type , time_type >( m_b2 , dt) );
    }

    template< class System , class StateIn , class DerivIn , class StateOut >
    void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt ,
            const time_type &t , StateOut &out , const time_type &dt )
    {
        typedef typename boost::unwrap_reference< System >::type unwrapped_system_type;
        unwrapped_system_type &sys = system;

        m_deriv_adjuster.adjust_size_by_policy( in , adjust_size_policy() );
        m_state_adjuster.adjust_size_by_policy( in , adjust_size_policy() );

        // actual calculation done in generic_rk.hpp
        m_rk_algorithm.do_step( sys , in , dxdt , t , out , dt , m_x_tmp , m_F );
    }

    template< class StateType >
    void adjust_size( const StateType &x )
    {
        m_deriv_adjuster.adjust_size( x );
        m_state_adjuster.adjust_size( x );
        stepper_base_type::adjust_size( x );
    }

private:

    rk_algorithm_type m_rk_algorithm;
    const coef_b_type m_b2;

    size_adjuster< deriv_type , StageCount-1 > m_deriv_adjuster;
    size_adjuster< state_type , 1 > m_state_adjuster;

    state_type m_x_tmp;
    deriv_type m_F[StageCount-1];

};

}
}
}


#endif /* EXPLICIT_ERROR_GENERIC_RK_HPP_ */
