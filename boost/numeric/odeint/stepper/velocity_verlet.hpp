/*
  [auto_generated]
  boost/numeric/odeint/stepper/velocity_verlet.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_VELOCITY_VERLET_HPP_DEFINED
#define BOOST_NUMERIC_ODEINT_STEPPER_VELOCITY_VERLET_HPP_DEFINED

#include <boost/numeric/odeint/stepper/base/algebra_stepper_base.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/operations_dispatcher.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>
#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>

#include <boost/numeric/odeint/util/bind.hpp>
#include <boost/numeric/odeint/util/copy.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>
// #include <boost/numeric/odeint/util/is_pair.hpp>
// #include <boost/array.hpp>



namespace boost {
namespace numeric {
namespace odeint {



/**
 * \class velocity_verlet
 * \brief The Velocity-Verlet algorithm.
 *
 * The Adams-Bashforth method is a multi-step algorithm with configurable step
 * number. The step number is specified as template parameter Steps and it 
 * then uses the result from the previous Steps steps. See also
 * <a href="http://en.wikipedia.org/wiki/Linear_multistep_method">en.wikipedia.org/wiki/Linear_multistep_method</a>.
 * Currently, a maximum of Steps=8 is supported.
 * The method is explicit and fulfills the Stepper concept. Step size control
 * or continuous output are not provided.
 * 
 * \tparam Steps The number of steps (maximal 8).
 * \tparam State The state type.
 * \tparam Value The value type.
 * \tparam Deriv The type representing the time derivative of the state.
 * \tparam Time The time representing the independent variable - the time.
 * \tparam Algebra The algebra type.
 * \tparam Operations The operations type.
 * \tparam Resizer The resizer policy type.
 * \tparam InitializingStepper The stepper for the first two steps.
 */

template <
    class Coor ,
    class Velocity = Coor ,
    class Value = double ,
    class Acceleration = Coor ,
    class Time = Value ,
    class TimeSq = Time ,
    class Algebra = typename algebra_dispatcher< Coor >::algebra_type ,
    class Operations = typename operations_dispatcher< Coor >::operations_type ,
    class Resizer = initially_resizer
      >
class velocity_verlet : public algebra_stepper_base< Algebra , Operations >
{
public:

    typedef algebra_stepper_base< Algebra , Operations > algebra_stepper_base_type;
    typedef typename algebra_stepper_base_type::algebra_type algebra_type;
    typedef typename algebra_stepper_base_type::operations_type operations_type;

    typedef Coor coor_type;
    typedef Velocity velocity_type;
    typedef Acceleration acceleration_type;
    typedef std::pair< coor_type , velocity_type > state_type;
    typedef std::pair< velocity_type , acceleration_type > deriv_type;
    typedef state_wrapper< acceleration_type > wrapped_acceleration_type;
    typedef Value value_type;
    typedef Time time_type;
    typedef TimeSq time_square_type;
    typedef Resizer resizer_type;
    typedef stepper_tag stepper_category;

    typedef unsigned short order_type;

    static const order_type order_value = 1;

    order_type order( void ) const
    {
        return order_value;
    }


    /**
     * \fn velocity_verlet::velocity_verlet( const algebra_type &algebra )
     * \brief Constructs the velocity_verlet class. This constructor can be used as a default
     * constructor if the algebra has a default constructor. 
     * \param algebra A copy of algebra is made and stored.
     */
    velocity_verlet( const algebra_type & algebra = algebra_type() )
        : algebra_stepper_base_type( algebra ) , m_first_call( true )
        , m_a1() , m_a2() , m_current_a1( true ) { }


    template< class System , class StateInOut >
    void do_step( System system , StateInOut & x , time_type t , time_type dt )
    {
        do_step_v1( system , x , t , dt );
    }
    
    template< class System , class StateInOut >
    void do_step( System system , const StateInOut & x , time_type t , time_type dt )
    {
        do_step_v1( system , x , t , dt );
    }

    template< class System , class CoorIn , class VelocityIn , class AccelerationIn ,
                             class CoorOut , class VelocityOut , class AccelerationOut >
    void do_step( System system , CoorIn const & qin , VelocityIn const & pin , AccelerationIn const & ain ,
                  CoorOut & qout , VelocityOut & pout , AccelerationOut & aout , time_type t , time_type dt )
    {
        const value_type one = static_cast< value_type >( 1.0 );
        const value_type one_half = static_cast< value_type >( 0.5 );

        algebra_stepper_base_type::m_algebra.for_each4(
            qout , qin , pin , ain ,
            typename operations_type::template scale_sum3< value_type , time_type , time_square_type >( one , one_half * dt , one * dt * dt ) );

        typename odeint::unwrap_reference< System >::type & sys = system;

        sys( qout , pin , aout , t + dt );

        algebra_stepper_base_type::m_algebra.for_each4(
            pout , pin , ain , aout ,
            typename operations_type::template scale_sum3< value_type , time_type , time_type >( one , one_half * dt , one_half * dt ) );
    }


    template< class StateIn >
    void adjust_size( const StateIn & x )
    {
        if( resize_impl( x ) )
            m_first_call = true;
    }

    void reset( void )
    {
        m_first_call = true;
    }

    template< class AccelerationIn >
    void initialize( const AccelerationIn & ain )
    {
        // alloc a
        m_resizer.adjust_size( ain ,
                               detail::bind( &velocity_verlet::template resize_impl< AccelerationIn > ,
                                             detail::ref( *this ) , detail::_1 ) );
        boost::numeric::odeint::copy( ain , get_current_acc() );
        m_first_call = false;
    }

    template< class System , class CoorIn , class VelocityIn >
    void initialize( System system , const CoorIn & qin , const VelocityIn & pin , time_type t )
    {
        m_resizer.adjust_size( qin ,
                               detail::bind( &velocity_verlet::template resize_impl< CoorIn > ,
                                             detail::ref( *this ) , detail::_1 ) );
        initialize_acc( system , qin , pin , t );
    }

    bool is_initialized( void ) const
    {
        return ! m_first_call;
    }


private:
    
    template< class System , class CoorIn , class VelocityIn >
    void initialize_acc( System system , const CoorIn & qin , const VelocityIn & pin , time_type t )
    {
        typename odeint::unwrap_reference< System >::type & sys = system;
        sys( qin , pin , get_current_acc() , t );
        m_first_call = false;
    }
    
    template< class System , class StateInOut >
    void do_step_v1( System system , StateInOut & x , time_type t , time_type dt )
    {
        typedef typename odeint::unwrap_reference< StateInOut >::type state_in_type;
        typedef typename odeint::unwrap_reference< typename state_in_type::first_type >::type coor_in_type;
        typedef typename odeint::unwrap_reference< typename state_in_type::second_type >::type momentum_in_type;
        
        typedef typename boost::remove_reference< coor_in_type >::type xyz_type;
        state_in_type & statein = x;
        coor_in_type & qinout = statein.first;
        momentum_in_type & pinout = statein.second;

        // alloc a
        if( m_resizer.adjust_size( qinout ,
                                   detail::bind( &velocity_verlet::template resize_impl< xyz_type > ,
                                                 detail::ref( *this ) , detail::_1 ) )
         || m_first_call )
        {
            initialize_acc( system , qinout , pinout , t );
        }

        // check first
        do_step( system , qinout , pinout , get_current_acc() , qinout , pinout , get_old_acc() , t , dt );
        toggle_current_acc();
    }

    template< class StateIn >
    bool resize_impl( const StateIn & x )
    {
        bool resized = false;
        resized |= adjust_size_by_resizeability( m_a1 , x , typename is_resizeable< acceleration_type >::type() );
        resized |= adjust_size_by_resizeability( m_a2 , x , typename is_resizeable< acceleration_type >::type() );
        return resized;
    }

    acceleration_type & get_current_acc( void )
    {
        return m_current_a1 ? m_a1.m_v : m_a2.m_v ;
    }

    const acceleration_type & get_current_acc( void ) const
    {
        return m_current_a1 ? m_a1.m_v : m_a2.m_v ;
    }

    acceleration_type & get_old_acc( void )
    {
        return m_current_a1 ? m_a2.m_v : m_a1.m_v ;
    }

    const acceleration_type & get_old_acc( void ) const
    {
        return m_current_a1 ? m_a2.m_v : m_a1.m_v ;
    }

    void toggle_current_acc( void )
    {
        m_current_a1 = ! m_current_a1;
    }

    resizer_type m_resizer;
    bool m_first_call;
    wrapped_acceleration_type m_a1 , m_a2;
    bool m_current_a1;
};




} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_VELOCITY_VERLET_HPP_DEFINED
