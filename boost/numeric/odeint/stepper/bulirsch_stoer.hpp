/*
 [auto_generated]
 boost/numeric/odeint/stepper/bulirsch_stoer.hpp

 [begin_description]
 Implementaiton of the Burlish-Stoer method.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_BULIRSCH_STOER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_BULIRSCH_STOER_HPP_INCLUDED


#include <iostream>

#include <algorithm>

#include <boost/ref.hpp>
#include <boost/bind.hpp>

#include <boost/numeric/odeint/stepper/controlled_error_stepper.hpp>
#include <boost/numeric/odeint/stepper/modified_midpoint.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
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
class bulirsch_stoer {

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
    typedef controlled_stepper_tag stepper_category;

    typedef bulirsch_stoer< State , Value , Deriv , Time , Algebra , Operations , Resizer > controlled_error_bs_type;

    typedef std::vector< time_type > value_vector;
    typedef std::vector< value_vector > value_matrix;
    typedef std::vector< size_t > int_vector;

    bulirsch_stoer(
            time_type eps_abs = 1E-6 , time_type eps_rel = 1E-6 ,
            time_type factor_x = 1.0 , time_type factor_dxdt = 1.0 )
    : m_error_checker( eps_abs , eps_rel , factor_x, factor_dxdt ),
      m_k_max(8) ,
      m_safety1(0.25) , m_safety2(0.7),
      m_max_dt_factor( 0.1 ) , m_min_step_scale(5E-5) , m_max_step_scale(0.7),
      m_last_step_rejected( false ) , m_first( true ) ,
      m_dt_last( 1.0E30 ) ,
      m_current_eps( -1.0 ) ,
      m_error( m_k_max ) ,
      m_a( m_k_max+1 ) ,
      m_alpha( m_k_max , value_vector( m_k_max ) ) ,
      m_interval_sequence( m_k_max+1 ) ,
      m_coeff( m_k_max+1 ) ,
      m_cost( m_k_max+1 ) ,
      m_times( m_k_max ) ,
      m_table( m_k_max ) ,
      STEPFAC1( 0.65 ) , STEPFAC2( 0.94 ) , STEPFAC3( 0.02 ) , STEPFAC4( 4.0 ) , KFAC1( 0.8 ) , KFAC2( 0.9 )
    {
        for( unsigned short i = 0; i < m_k_max+1; i++ )
        {
            m_interval_sequence[i] = 2 * (i+1);
            if( i == 0 )
                m_cost[i] = m_interval_sequence[i];
            else
                m_cost[i] = m_cost[i-1] + m_interval_sequence[i];
            m_coeff[i].resize(i);
            for( size_t k = 0 ; k < i ; ++k  )
            {
                const time_type r = static_cast< time_type >( m_interval_sequence[i] ) / static_cast< time_type >( m_interval_sequence[k] );
                m_coeff[i][k] = 1.0 / ( r*r - static_cast< time_type >( 1.0 ) ); // coefficients for extrapolation
                //std::cout << i << "," << k << " " << m_coeff[i][k] << '\t' ;
            }
            //std ::cout << std::endl;
            // crude estimate of optimal order
            const time_type logfact( -log10( std::max( eps_rel , 1.0E-12 ) ) * 0.6 + 0.5 );
            m_current_k_opt = std::max( 1 , std::min( static_cast<int>( m_k_max-1 ) , static_cast<int>( logfact ) ));
            //m_current_k_opt = m_k_max - 1;
            //std::cout << m_cost[i] << std::endl;
        }

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
        static const time_type val1( static_cast< time_type >( 1.0 ) );

        typename boost::unwrap_reference< System >::type &sys = system;
        if( m_resizer.adjust_size( in , boost::bind( &controlled_error_bs_type::resize< StateIn > , boost::ref( *this ) , _1 ) ) )
            reset(); // system resized -> reset
        if( dt != m_dt_last )
            reset(); // step size changed from outside -> reset

        bool reject( true );
        m_dt_last = dt;

        value_vector h_opt( m_k_max+1 );
        value_vector work( m_k_max+1 );

        size_t k_final = 0;

        //std::cout << "t=" << t <<", dt=" << dt << ", k_opt=" << m_current_k_opt << std::endl;

        for( size_t k = 0 ; k <= m_current_k_opt+1 ; k++ )
        {
            //std::cout << "k=" << k <<": " << std::endl;
            m_midpoint.set_steps( m_interval_sequence[k] );
            if( k == 0 )
            {
                m_midpoint.do_step( sys , in , dxdt , t , out , dt );
            }
            else
            {
                m_midpoint.do_step( sys , in , dxdt , t , m_table[k-1].m_v , dt );
                extrapolate( k , out );
                // get error estimate
                m_algebra.for_each3( m_err.m_v , out , m_table[0].m_v ,
                        typename operations_type::template scale_sum2< time_type , time_type >( val1 , -val1 ) );
                const time_type error = m_error_checker.error( m_algebra , in , dxdt , m_err.m_v , dt );
                h_opt[k] = calc_h_opt( dt , error , k );
                work[k] = m_cost[k]/h_opt[k];
                //std::cout << '\t' << "h_opt=" << h_opt[k] << ", work=" << work[k] << std::endl;
                //std::cout << '\t' << "error: " << error << std::endl;
                if( m_first && (error <= static_cast< time_type >( 1.0 )) )
                { // this is the first step, convergence does not have to be in order window
                    //std::cout << '\t' << "convergence in first step" << std::endl;
                    reject = false;
                    k_final = k;
                    break; // leave k-loop
                }
                if( in_convergence_window( k ) )
                {
                    k_final = k;
                    if( error <= 1.0 ) // convergence in convergence window
                    {
                        reject = false;
                        break; // we are good, leave k-loop
                    }
                    else
                    {
                        if( should_reject( error , k ) )
                        {
                            reject = true;
                            break;
                        }
                    }
                }
            }
        }

        //calculate optimal order and stepsize
        controlled_step_result result = set_k_opt( k_final , work , h_opt , dt );

        m_last_step_rejected = reject;
        if( reject )
            result = step_size_decreased;
        else
            t += m_dt_last;

        m_dt_last = dt;

        m_first = false;

        return result;
    }

    void reset()
    {
        m_first = true;
        m_last_step_rejected = false;
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
    bool resize( const StateIn &x )
    {
        bool resized( false );
        for( size_t i = 0 ; i < m_k_max ; ++i )
            resized |= adjust_size_by_resizeability( m_table[i] , x , typename wrapped_state_type::is_resizeable() );
        resized |= adjust_size_by_resizeability( m_err , x , typename wrapped_state_type::is_resizeable() );
        return resized;
    }

    template< class StateIn >
    void adjust_size( const StateIn &x )
    {
        resize_m_dxdt( x );
        resize_m_xnew( x );
        resize( x );
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


    template< class StateInOut >
    void extrapolate( const size_t k , StateInOut &xest )
    //polynomial extrapolation, see http://www.nr.com/webnotes/nr3web21.pdf
    {
        //std::cout << "extrapolate k=" << k << ":" << std::endl;
        static const time_type val1 = static_cast< time_type >( 1.0 );
        for( int j=k-1 ; j>0 ; --j )
        {
            //std::cout << '\t' << m_coeff[k][j];
            m_algebra.for_each3( m_table[j-1].m_v , m_table[j].m_v , m_table[j-1].m_v ,
                    typename operations_type::template scale_sum2< time_type , time_type >( val1 + m_coeff[k][j] , -m_coeff[k][j] ) );
        }
        //std::cout << std::endl << m_coeff[k][0] << std::endl;
        m_algebra.for_each3( xest , m_table[0].m_v , xest ,
                typename operations_type::template scale_sum2< time_type , time_type >( val1 + m_coeff[k][0] , -m_coeff[k][0]) );
    }

    time_type calc_h_opt( const time_type h , const value_type error , const size_t k ) const
    {
        time_type expo=1.0/(2*k+1);
        time_type facmin = std::pow( STEPFAC3 , expo );
        time_type fac;
        if (error == 0.0)
            fac=1.0/facmin;
        else
        {
            fac = STEPFAC2 / std::pow( error / STEPFAC1 , expo );
            fac = std::max( facmin/STEPFAC4 , std::min( 1.0/facmin , fac ) );
        }
        return std::abs(h*fac);
    }

    controlled_step_result set_k_opt( const size_t k , const value_vector &work , const value_vector &h_opt , time_type &dt )
    {
        //std::cout << "finding k_opt..." << std::endl;
        if( k == 1 )
        {
            m_current_k_opt = 2;
            //dt = h_opt[ m_current_k_opt-1 ] * m_cost[ m_current_k_opt ] / m_cost[ m_current_k_opt-1 ] ;
            return success_step_size_increased;
        }
        if( (work[k-1] < KFAC1*work[k]) || (k == m_k_max) )
        {   // order decrease
            m_current_k_opt = k-1;
            dt = h_opt[ m_current_k_opt ];
            return success_step_size_increased;
        }
        else if( (work[k] < KFAC2*work[k-1]) || m_last_step_rejected || (k == m_k_max-1) )
        {   // same order - also do this if last step got rejected
            m_current_k_opt = k;
            dt = h_opt[ m_current_k_opt ];
            return success_step_size_unchanged;
        }
        else
        {   // order increase - only if last step was not rejected
            m_current_k_opt = k+1;
            dt = h_opt[ m_current_k_opt-1 ] * m_cost[ m_current_k_opt ] / m_cost[ m_current_k_opt-1 ] ;
            return success_step_size_increased;
        }
    }

    bool in_convergence_window( const size_t k ) const
    {
        if( (k == m_current_k_opt-1) && !m_last_step_rejected )
            return true; // decrease stepsize only if last step was not rejected
        return ( (k == m_current_k_opt) || (k == m_current_k_opt+1) );
    }

    bool should_reject( const time_type error , const size_t k ) const
    {
        if( (k == m_current_k_opt-1) )
        {
            const time_type d = m_interval_sequence[m_current_k_opt] * m_interval_sequence[m_current_k_opt+1] /
                    (m_interval_sequence[0]*m_interval_sequence[0]);
            //step will fail, criterion 17.3.17 in NR
            return ( error > d*d );
        }
        else if( k == m_current_k_opt )
        {
            const time_type d = m_interval_sequence[m_current_k_opt] / m_interval_sequence[0];
            return ( error > d*d );
        } else
            return error > 1.0;
    }

    default_error_checker< value_type, algebra_type , operations_type > m_error_checker;
    modified_midpoint< state_type , value_type , deriv_type , time_type , algebra_type , operations_type , resizer_type > m_midpoint;

    const size_t m_k_max;

    const time_type m_safety1;
    const time_type m_safety2;
    const time_type m_max_dt_factor;
    const time_type m_min_step_scale;
    const time_type m_max_step_scale;

    bool m_last_step_rejected;
    bool m_first;

    time_type m_dt_last;
    time_type m_t_last;
    time_type m_current_eps;

    size_t m_current_k_max;
    size_t m_current_k_opt;

    algebra_type m_algebra;

    resizer_type m_dxdt_resizer;
    resizer_type m_xnew_resizer;
    resizer_type m_resizer;

    wrapped_state_type m_xnew;
    wrapped_state_type m_err;
    wrapped_deriv_type m_dxdt;

    value_vector m_error; // errors of repeated midpoint steps and extrapolations
    value_vector m_a; // stores the work (number of f calls) required for the orders
    value_matrix m_alpha; // stores convergence factor for stepsize adjustment
    int_vector m_interval_sequence; // stores the successive interval counts
    value_matrix m_coeff;
    int_vector m_cost; // costs for interval count

    value_vector m_times;
    std::vector< wrapped_state_type > m_table; // sequence of states for extrapolation

    const time_type STEPFAC1 , STEPFAC2 , STEPFAC3 , STEPFAC4 , KFAC1 , KFAC2;
};

}
}
}

#endif // BOOST_NUMERIC_ODEINT_STEPPER_BULIRSCH_STOER_HPP_INCLUDED
