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

template< class Value = double >
struct rk4_coefficients_a1 : boost::array< Value , 1 >
{
    rk4_coefficients_a1( void )
    {
        (*this)[0] = static_cast< Value >( 0.5 );
    }
};

template< class Value = double >
struct rk4_coefficients_a2 : boost::array< Value , 2 >
{
    rk4_coefficients_a2( void )
    {
        (*this)[0] = static_cast<Value>(0.0);
        (*this)[1] = static_cast<Value>(0.5);
    }
};


template< class Value = double >
struct rk4_coefficients_a3 : boost::array< Value , 3 >
{
    rk4_coefficients_a3( void )
    {
        (*this)[0] = static_cast<Value>(0.0);
        (*this)[1] = static_cast<Value>(0.0);
        (*this)[2] = static_cast<Value>(1.0);
    }
};

template< class Value = double >
struct rk4_coefficients_b : boost::array< Value , 4 >
{
    rk4_coefficients_b( void )
    {
        (*this)[0] = static_cast<Value>(1.0)/static_cast<Value>(6.0);
        (*this)[1] = static_cast<Value>(1.0)/static_cast<Value>(3.0);
        (*this)[2] = static_cast<Value>(1.0)/static_cast<Value>(3.0);
        (*this)[3] = static_cast<Value>(1.0)/static_cast<Value>(6.0);
    }
};

template< class Value = double >
struct rk4_coefficients_c : boost::array< Value , 4 >
{
    rk4_coefficients_c( void )
    {
        (*this)[0] = static_cast<Value>(0.0);
        (*this)[1] = static_cast<Value>(0.5);
        (*this)[2] = static_cast<Value>(0.5);
        (*this)[3] = static_cast<Value>(1.0);
    }
};


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
            fusion::make_vector( rk4_coefficients_a1<Value>() , rk4_coefficients_a2<Value>() , rk4_coefficients_a3<Value>() ) ,
            rk4_coefficients_b<Value>() , rk4_coefficients_c<Value>() )
    { }

};

}
}
}


#endif /* EXPLICIT_RK4_GENERIC_HPP_ */
