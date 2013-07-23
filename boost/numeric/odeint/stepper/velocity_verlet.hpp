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

// #include <boost/numeric/odeint/util/bind.hpp>
// #include <boost/numeric/odeint/util/copy.hpp>
// #include <boost/numeric/odeint/util/is_pair.hpp>
// #include <boost/numeric/odeint/util/resizer.hpp>
// #include <boost/array.hpp>


namespace boost {
namespace numeric {
namespace odeint {

template<
    class Coor ,
    class Momentum = Coor ,
    class Value = double ,
    class CoorDeriv = Coor ,
    class MomentumDeriv = Coor ,
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
    typedef Momentum momentum_type;
    typedef std::pair< coor_type , momentum_type > state_type;
    typedef CoorDeriv coor_deriv_type;
    typedef MomentumDeriv momentum_deriv_type;
    typedef std::pair< coor_deriv_type , momentum_deriv_type > deriv_type;    
    typedef state_wrapper< coor_deriv_type> wrapped_coor_deriv_type;
    typedef state_wrapper< momentum_deriv_type > wrapped_momentum_deriv_type;
    typedef Value value_type;
    typedef Time time_type;
    typedef TimeSq time_square_type;
    typedef Resizer resizer_type;
    typedef stepper_tag stepper_category;
    
    typedef unsigned short order_type;

    static const order_type order_value = 1;
    
    
    template< class System , class StateIn >
    void do_step( System system , StateIn &x , time_type t , time_type dt )
    {
    }
    
    template< class System , class CoorIn , class MomentumIn , class MomentumDerivIn , class CoorOut , class MomentumOut , class MomentumDerivOut >
    void do_step( System system , CoorIn const& qin , MomentumIn const& pin , MomentumDerivIn const& ain , CoorOut &qout , MomentumOut &pout , MomentumDerivOut &aout , time_type t , time_type dt )
    {
        const value_type one = static_cast< value_type >( 1.0 );
        const value_type one_half = static_cast< value_type >( 0.5 );
        
        algebra_stepper_base_type::for_each4( qout , qin , pin , ain , 
                                              typename operations_type::template scale_sum3< value_type , time_type , time_square_type >( one , one_half * dt , one * dt * dt ) );
        
        typename odeint::unwrap_reference< System >::type &sys = system;
        
        sys( qout , pin , aout );
        
        algebra_stepper_base_type::for_each4( pout , pin , ain , aout ,
                                              typename operations_type::template scale_sum3< value_type , time_type , time_type >( one , one_half * dt , one_half * dt ) );
    }
    
private:
    
    wrapped_momentum_deriv_type m_dqdt;
};




} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_VELOCITY_VERLET_HPP_DEFINED
