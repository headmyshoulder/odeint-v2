/*
 * explicit_error_rk54_ck_generic.hpp
 *
 *  Created on: Jun 3, 2011
 *      Author: mario
 */

#include <boost/fusion/container.hpp>

#include <boost/numeric/odeint/stepper/explicit_error_generic_rk.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>

#include <boost/array.hpp>


namespace fusion = boost::fusion;


namespace boost {
namespace numeric {
namespace odeint {

template< class Value = double >
struct rk54_ck_coefficients_a1 : boost::array< Value , 1 >
{
    rk54_ck_coefficients_a1( void )
    {
        (*this)[0] = static_cast< Value >( 1 )/static_cast< Value >( 5 );
    }
};

template< class Value = double >
struct rk54_ck_coefficients_a2 : boost::array< Value , 2 >
{
    rk54_ck_coefficients_a2( void )
    {
        (*this)[0] = static_cast<Value>( 3 )/static_cast<Value>( 40 );
        (*this)[1] = static_cast<Value>( 9 )/static_cast<Value>( 40 );
    }
};


template< class Value = double >
struct rk54_ck_coefficients_a3 : boost::array< Value , 3 >
{
    rk54_ck_coefficients_a3( void )
    {
        (*this)[0] = static_cast<Value>( 3 )/static_cast<Value>( 10 );
        (*this)[1] = static_cast<Value>( -9 )/static_cast<Value>( 10 );
        (*this)[2] = static_cast<Value>( 6 )/static_cast<Value>( 5 );
    }
};

template< class Value = double >
struct rk54_ck_coefficients_a4 : boost::array< Value , 4 >
{
    rk54_ck_coefficients_a4( void )
    {
        (*this)[0] = static_cast<Value>( -11 )/static_cast<Value>( 54 );
        (*this)[1] = static_cast<Value>( 5 )/static_cast<Value>( 2 );
        (*this)[2] = static_cast<Value>( -70 )/static_cast<Value>( 27 );
        (*this)[3] = static_cast<Value>( 35 )/static_cast<Value>( 27 );
    }
};

template< class Value = double >
struct rk54_ck_coefficients_a5 : boost::array< Value , 5 >
{
    rk54_ck_coefficients_a5( void )
    {
        (*this)[0] = static_cast<Value>( 1631 )/static_cast<Value>( 55296 );
        (*this)[1] = static_cast<Value>( 175 )/static_cast<Value>( 512 );
        (*this)[2] = static_cast<Value>( 575 )/static_cast<Value>( 13824 );
        (*this)[3] = static_cast<Value>( 44275 )/static_cast<Value>( 110592 );
        (*this)[4] = static_cast<Value>( 253 )/static_cast<Value>( 4096 );
    }
};

template< class Value = double >
struct rk54_ck_coefficients_b : boost::array< Value , 6 >
{
    rk54_ck_coefficients_b( void )
    {
        (*this)[0] = static_cast<Value>( 37 )/static_cast<Value>( 378 );
        (*this)[1] = static_cast<Value>( 0 );
        (*this)[2] = static_cast<Value>( 250 )/static_cast<Value>( 621 );
        (*this)[3] = static_cast<Value>( 125 )/static_cast<Value>( 594 );
        (*this)[4] = static_cast<Value>( 0 );
        (*this)[5] = static_cast<Value>( 512 )/static_cast<Value>( 1771 );
    }
};

template< class Value = double >
struct rk54_ck_coefficients_db : boost::array< Value , 6 >
{
    rk54_ck_coefficients_db( void )
    {
        (*this)[0] = static_cast<Value>( 37 )/static_cast<Value>( 378 ) - static_cast<Value>( 2825 )/static_cast<Value>( 27648 );
        (*this)[1] = static_cast<Value>( 0 );
        (*this)[2] = static_cast<Value>( 250 )/static_cast<Value>( 621 ) - static_cast<Value>( 18575 )/static_cast<Value>( 48384 );
        (*this)[3] = static_cast<Value>( 125 )/static_cast<Value>( 594 ) - static_cast<Value>( 13525 )/static_cast<Value>( 55296 );
        (*this)[4] = static_cast<Value>( -277 )/static_cast<Value>( 14336 );
        (*this)[5] = static_cast<Value>( 512 )/static_cast<Value>( 1771 ) - static_cast<Value>( 1 )/static_cast<Value>( 4 );
    }
};


template< class Value = double >
struct rk54_ck_coefficients_c : boost::array< Value , 6 >
{
    rk54_ck_coefficients_c( void )
    {
        (*this)[0] = static_cast<Value>(0);
        (*this)[1] = static_cast<Value>( 1 )/static_cast<Value>( 5 );
        (*this)[2] = static_cast<Value>( 3 )/static_cast<Value>( 10 );
        (*this)[3] = static_cast<Value>( 3 )/static_cast<Value>( 5 );
        (*this)[4] = static_cast<Value>( 3 )/static_cast<Value>( 5 );
        (*this)[4] = static_cast<Value>( 1 );
        (*this)[4] = static_cast<Value>( 7 )/static_cast<Value>( 8 );
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
class explicit_error_rk54_ck_generic : public explicit_error_generic_rk< 6 , 5 , 5 , 4 ,
        State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy >
{

public:

    typedef explicit_error_generic_rk< 6 , 5 , 5 , 4 , State , Value , Deriv , Value ,
                               Algebra , Operations , AdjustSizePolicy > stepper_base_type;

    typedef typename stepper_base_type::state_type state_type;
    typedef typename stepper_base_type::value_type value_type;
    typedef typename stepper_base_type::deriv_type deriv_type;
    typedef typename stepper_base_type::time_type time_type;
    typedef typename stepper_base_type::algebra_type algebra_type;
    typedef typename stepper_base_type::operations_type operations_type;
    typedef typename stepper_base_type::adjust_size_policy adjust_size_policy;
    typedef typename stepper_base_type::stepper_type stepper_type;

    explicit_error_rk54_ck_generic( void ) : stepper_base_type(
            fusion::make_vector( rk54_ck_coefficients_a1<Value>() ,
                                 rk54_ck_coefficients_a2<Value>() ,
                                 rk54_ck_coefficients_a3<Value>() ,
                                 rk54_ck_coefficients_a4<Value>() ,
                                 rk54_ck_coefficients_a5<Value>() ) ,
            rk54_ck_coefficients_b<Value>() , rk54_ck_coefficients_db<Value>() , rk54_ck_coefficients_c<Value>() )
    { }
};

}
}
}
