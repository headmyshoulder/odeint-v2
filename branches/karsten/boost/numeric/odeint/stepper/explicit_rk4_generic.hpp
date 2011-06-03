/*
 * explicit_rk4_generic.hpp
 *
 *  Created on: Jun 1, 2011
 *      Author: mario
 */

#ifndef EXPLICIT_RK4_GENERIC_HPP_
#define EXPLICIT_RK4_GENERIC_HPP_

#include <boost/fusion/container.hpp>

#include <boost/numeric/odeint/stepper/explicit_generic_rk.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>

#include <boost/array.hpp>


namespace fusion = boost::fusion;


namespace boost {
namespace numeric {
namespace odeint {

namespace constants_generic {

    const boost::array< double , 1 > rk4_a1 = {{ 0.5 }};
    const boost::array< double , 2 > rk4_a2 = {{ 0.0 , 0.5 }};
    const boost::array< double , 3 > rk4_a3 = {{ 0.0 , 0.0 , 1.0 }};

    const boost::array< double , 4 > rk4_b = {{ 1.0/6 , 1.0/3 , 1.0/3 , 1.0/6 }};
    const boost::array< double , 4 > rk4_c = {{ 0.0 , 0.5 , 0.5 , 1.0 }};

}

template<
    class State ,
    class Value = double ,
    class Deriv = State ,
    class Time = Value ,
    class Algebra = range_algebra ,
    class Operations = default_operations ,
    class AdjustSizePolicy = adjust_size_initially_tag
    >
class explicit_rk4_generic : public explicit_generic_rk< 4 , 4 , State , Value , Deriv , Value ,
                                                          Algebra , Operations , AdjustSizePolicy >
{

public:

    typedef explicit_generic_rk< 4 , 4 , State , Value , Deriv , Value ,
                               Algebra , Operations , AdjustSizePolicy > stepper_base_type;

    typedef typename stepper_base_type::state_type state_type;
    typedef typename stepper_base_type::value_type value_type;
    typedef typename stepper_base_type::deriv_type deriv_type;
    typedef typename stepper_base_type::time_type time_type;
    typedef typename stepper_base_type::algebra_type algebra_type;
    typedef typename stepper_base_type::operations_type operations_type;
    typedef typename stepper_base_type::adjust_size_policy adjust_size_policy;
    typedef typename stepper_base_type::stepper_type stepper_type;

    explicit_rk4_generic( void ) : stepper_base_type(
            fusion::make_vector( constants_generic::rk4_a1 , constants_generic::rk4_a2 , constants_generic::rk4_a3 ) ,
            constants_generic::rk4_b , constants_generic::rk4_c )
    { }

};

}
}
}


#endif /* EXPLICIT_RK4_GENERIC_HPP_ */
