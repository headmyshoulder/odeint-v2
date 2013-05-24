/*
  [auto_generated]
  boost/numeric/odeint/stepper/parallel_adams_bashforth_stepper.hpp

  [begin_description]
  parallel adams bashforth stepper, see
  "Execution Schemes for Parallel Adams Methods"
  Thomas Rauber and Gudula Runger
  [end_description]

  Copyright 2009-2013 Karsten Ahnert
  Copyright 2009-2013 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_PARALLEL_ADAMS_BASHFORTH_STEPPER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_PARALLEL_ADAMS_BASHFORTH_STEPPER_HPP_INCLUDED


#include <iostream>

#include <algorithm>

#include <boost/config.hpp> // for min/max guidelines
#include <boost/static_assert.hpp>

#include <boost/array.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/odeint/util/bind.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>

#include <boost/numeric/odeint/stepper/base/explicit_stepper_base.hpp>
#include <boost/numeric/odeint/stepper/parallel_extrapolation_stepper.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/operations_dispatcher.hpp>
#include <boost/numeric/odeint/stepper/detail/generic_rk_call_algebra.hpp>
#include <boost/numeric/odeint/stepper/detail/generic_rk_operations.hpp>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>
#include <boost/numeric/odeint/util/unit_helper.hpp>
#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>

#include <boost/numeric/odeint/stepper/detail/parallel_adams_bashforth_coefficients.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template<
    unsigned short Stages ,
    class State ,
    class Value = double ,
    class Deriv = State ,
    class Time = Value ,
    class Algebra = typename algebra_dispatcher< State >::algebra_type ,
    class Operations = typename operations_dispatcher< State >::operations_type ,
    class Resizer = initially_resizer
    >
#ifndef DOXYGEN_SKIP
class parallel_adams_bashforth_stepper : public explicit_stepper_base<
    parallel_adams_bashforth_stepper< Stages , State , Value , Deriv , Time , Algebra , Operations , Resizer > ,
    Stages+1 , State , Value , Deriv , Time , Algebra , Operations , Resizer >
#else
class parallel_adams_bashforth_stepper : public explicit_stepper_base
#endif
{

 private:
    // 2 < Stages < 6
    BOOST_STATIC_ASSERT( Stages > 2 );
    BOOST_STATIC_ASSERT( Stages < 6 );

 public:

#ifndef DOXYGEN_SKIP
    typedef explicit_stepper_base< parallel_adams_bashforth_stepper< Stages , State , Value , Deriv , Time , Algebra , Operations , Resizer > ,
            Stages+1 , State , Value , Deriv , Time , Algebra , Operations , Resizer > stepper_base_type;
#else
    typedef explicit_stepper_base< parallel_adams_bashforth_stepper< ... > , ... > stepper_base_type;
#endif

    typedef typename stepper_base_type::state_type state_type;
    typedef typename stepper_base_type::value_type value_type;
    typedef typename stepper_base_type::deriv_type deriv_type;
    typedef typename stepper_base_type::time_type time_type;
    typedef typename stepper_base_type::algebra_type algebra_type;
    typedef typename stepper_base_type::operations_type operations_type;
    typedef typename stepper_base_type::resizer_type resizer_type;

#ifndef DOXYGEN_SKIP
    typedef typename stepper_base_type::stepper_type stepper_type;
    typedef typename stepper_base_type::wrapped_state_type wrapped_state_type;
    typedef typename stepper_base_type::wrapped_deriv_type wrapped_deriv_type;

    typedef boost::numeric::ublas::c_matrix< value_type , Stages , Stages > value_matrix;

    // use extrapolation stepper with even order >= Stages+1 as initializing stepper
    typedef parallel_extrapolation_stepper< 2*( (Stages+2)/2 ) , State , Value , Deriv , Time , Algebra , Operations , Resizer > init_stepper_type;

#endif //DOXYGEN_SKIP

    typedef unsigned short order_type;
    static const order_type order_value = stepper_base_type::order_value;
    
    parallel_adams_bashforth_stepper( const algebra_type &algebra = algebra_type() )
        : stepper_base_type( algebra ) , m_init_state_count(0)
    {
        // initialize steppers
        for( size_t n=0 ; n<Stages-1 ; ++n )
        {
            m_init_steppers[n] = init_stepper_type( algebra );
        }
        // calculate coefficients from lobatto points
        detail::parallel_adams_bashforth_coefficients< value_type , Stages > a;
        value_matrix v_a , w_b , s;
        for( size_t n=0 ; n<Stages ; ++n )
        {
            v_a(n,0) = a[n];
            w_b(n,0) = 1.0;
            value_type tmp = 1.0;
            for( size_t j=1 ; j<Stages ; ++j )
            {
                v_a(n,j) = v_a(n,j-1) * (a[n]);
                tmp *= a[n]-1.0;
                w_b(n,j) = (j+1) * tmp;
            }
        }

        using boost::numeric::ublas::prod;
        value_matrix w_b_inv;
        invert_matrix( w_b , w_b_inv );
        s = prod( v_a , w_b_inv );


        for( size_t n=0 ; n<Stages ; ++n )
        {
            // generic_rk_scale_sum requires this coefficient layout
            m_coeff[n][0] = s(n,Stages-1);
            if( n<Stages-1 )
                m_init_coeff[n] = s(n,0) - 1.0;
            for( size_t j=1 ; j<Stages ; ++j )
            {
                m_coeff[n][j] = s(n,j-1);
                if( n<Stages-1 )
                    m_init_coeff[n] += s(n,j);
            }
        }
    }
    
    
    template< class System , class StateIn , class DerivIn , class StateOut >
    void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt ,
                       time_type t , StateOut &out , time_type dt )
    {
        m_resizer.adjust_size( in , detail::bind( &stepper_type::template resize_impl<StateIn> , detail::ref( *this ) , detail::_1 ) );

        while( m_init_state_count < Stages-1 )
        {
            m_init_steppers[m_init_state_count].do_step( system , in , dxdt , t ,
                    m_states[m_init_state_count].m_v , m_init_coeff[m_init_state_count]*dt );
            m_init_state_count++;
        }

        typename odeint::unwrap_reference< System >::type &sys = system;

        for( size_t n=0 ; n<Stages-1 ; ++n )
        {
            // compute derivative
            sys( m_states[n].m_v , m_derivs[n].m_v , t );
        }
        for( size_t n=0 ; n<Stages-1 ; ++n )
        {
            // calculates states
            detail::template generic_rk_call_algebra< Stages , Algebra >()(
                    this->m_algebra , m_states[n].m_v , in , dxdt , m_derivs ,
                    detail::generic_rk_scale_sum< Stages , Operations , Value , Time >( m_coeff[n] , dt) );
        }
        // calculate last line
        detail::template generic_rk_call_algebra< Stages , Algebra >()(
                this->m_algebra , out , in , dxdt , m_derivs ,
                detail::generic_rk_scale_sum< Stages , Operations , Value , Time >( m_coeff[Stages-1] , dt) );
    }

    template< class StateIn >
    void adjust_size( const StateIn &x )
    {
        resize_impl( x );
    }

private:

    template< class StateIn >
    bool resize_impl( const StateIn &x )
    {
        bool resized( false );
        for( size_t i = 0 ; i < Stages-1 ; ++i )
        {
            resized |= adjust_size_by_resizeability( m_states[i] , x , typename is_resizeable<state_type>::type() );
            resized |= adjust_size_by_resizeability( m_derivs[i] , x , typename is_resizeable<state_type>::type() );
        }
        resized |= adjust_size_by_resizeability( m_deriv , x , typename is_resizeable<state_type>::type() );
        return resized;
    }

    void invert_matrix( value_matrix &m , value_matrix &inv )
    {
        using namespace boost::numeric::ublas;
        typedef permutation_matrix<std::size_t> pmatrix;
        // create a permutation matrix for the LU-factorization
        pmatrix pm(m.size1());
        // perform LU-factorization
        int res = lu_factorize(m,pm);
            if( res != 0 ) assert("matric inversion failed in parallel_adams_bashforht, shouldn't happen?!");
        // create identity matrix of "inverse"
        inv.assign( identity_matrix<value_type>(m.size1()) );
        // backsubstitute to get the inverse
        lu_substitute(m, pm, inv);
    }

 public:
    wrapped_state_type m_states[Stages-1];

 private:
    resizer_type m_resizer;
    
    wrapped_deriv_type m_derivs[Stages-1];
    wrapped_deriv_type m_deriv;

    init_stepper_type m_init_steppers[Stages-1];

    boost::array<value_type,Stages> m_coeff[Stages];
    value_type m_init_coeff[Stages-1];

    int m_init_state_count;
};

/******** DOXYGEN *******/

/**
 * \class extrapolation_stepper
 * \brief Extrapolation stepper with configurable order, and error estimation.
 *
 * The extrapolation stepper is a stepper with error estimation and configurable order. The order is given as
 * template parameter and needs to be an _odd_ number. The stepper is based on several executions of the
 * modified midpoint method and a Richardson extrapolation. This is essentially the same technique as for
 * bulirsch_stoer, but without the variable order.
 *
 * \note The Order parameter has to be an even number.
 */

} } }
#endif
