/*
 [auto_generated]
 boost/numeric/odeint/stepper/symplectic_rkn_sb3a_mclachlan.hpp
 
 [begin_description]
 Implementation of the symplectic MacLachlan stepper for separable Hamiltonian system.
 [end_description]
 
 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky
 
 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_SYMPLECTIC_RKN_SB3A_MCLACHLAN_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_SYMPLECTIC_RKN_SB3A_MCLACHLAN_HPP_INCLUDED


#include <boost/numeric/odeint/stepper/base/symplectic_rkn_stepper_base.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>

#include <boost/array.hpp>

namespace boost {
namespace numeric {
namespace odeint {

namespace detail {
namespace symplectic_rkn_sb3a_mclachlan {

    /*
      exp( a1 t A ) exp( b1 t B )
      exp( a2 t A ) exp( b2 t B )
      exp( a3 t A ) exp( b3 t B ) exp( a3 t A )
      exp( b2 t B ) exp( a2 t A )
      exp( b1 t B ) exp( a1 t A )
    */

    template< class Value >
    struct coef_a_type : public boost::array< Value , 6 >
    {
        coef_a_type( void )
        {
            (*this)[0] = static_cast< Value >( 0.40518861839525227722 );
            (*this)[1] = static_cast< Value >( -0.28714404081652408900 );
            (*this)[2] = static_cast< Value >( 1 ) / static_cast< Value >( 2 ) - ( (*this)[0] + (*this)[1] );
            (*this)[3] = (*this)[2];
            (*this)[4] = (*this)[1];
            (*this)[5] = (*this)[0];

        }
    };

    template< class Value >
    struct coef_b_type : public boost::array< Value , 6 >
    {
        coef_b_type( void )
        {
            (*this)[0] = static_cast< Value >( -3 ) / static_cast< Value >( 73 );
            (*this)[1] = static_cast< Value >( 17 ) / static_cast< Value >( 59 );
            (*this)[2] = static_cast< Value >( 1 ) - static_cast< Value >( 2 ) * ( (*this)[0] + (*this)[1] );
            (*this)[3] = (*this)[1];
            (*this)[4] = (*this)[0];
            (*this)[5] = static_cast< Value >( 0 );
        }
    };

} // namespace symplectic_rkn_sb3a_mclachlan
} // namespace detail



template<
    class Coor ,
    class Momentum = Coor ,
    class Value = double ,
    class CoorDeriv = Coor ,
    class MomentumDeriv = Coor ,
    class Time = Value ,
    class Algebra = range_algebra ,
    class Operations = default_operations ,
    class Resizer = initially_resizer
    >
class symplectic_rkn_sb3a_mclachlan :
        public symplectic_nystroem_stepper_base
<
    6 , 4 ,
    symplectic_rkn_sb3a_mclachlan< Coor , Momentum , Value , CoorDeriv , MomentumDeriv , Time , Algebra , Operations , Resizer > ,
    Coor , Momentum , Value , CoorDeriv , MomentumDeriv , Time , Algebra , Operations , Resizer
    >
{
public:

    typedef symplectic_nystroem_stepper_base
    <
    6 , 4 ,
    symplectic_rkn_sb3a_mclachlan< Coor , Momentum , Value , CoorDeriv , MomentumDeriv , Time , Algebra , Operations , Resizer > ,
    Coor , Momentum , Value , CoorDeriv , MomentumDeriv , Time , Algebra , Operations , Resizer
    > stepper_base_type;

    typedef typename stepper_base_type::algebra_type algebra_type;
    typedef typename stepper_base_type::value_type value_type;

    symplectic_rkn_sb3a_mclachlan( const algebra_type &algebra = algebra_type() )
        : stepper_base_type(
            detail::symplectic_rkn_sb3a_mclachlan::coef_a_type< value_type >() ,
            detail::symplectic_rkn_sb3a_mclachlan::coef_b_type< value_type >() ,
            algebra )
    { }
};



} // namespace odeint
} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_ODEINT_STEPPER_SYMPLECTIC_RKN_SB3A_MCLACHLAN_HPP_INCLUDED
