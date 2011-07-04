/*
 * symplectic_euler.hpp
 *
 *  Created on: Feb 12, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_SYMPLECTIC_EULER_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_SYMPLECTIC_EULER_HPP_

#include <boost/numeric/odeint/stepper/base/symplectic_rkn_stepper_base.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>

#include <boost/numeric/odeint/stepper/detail/macros.hpp>

#include <boost/array.hpp>

namespace boost {
namespace numeric {
namespace odeint {

namespace detail {
namespace symplectic_euler_coef {

	template< class Value >
	struct coef_a_type : public boost::array< Value , 1 >
	{
		coef_a_type( void )
		{
			(*this)[0] = 1.0;
		}
	};

	template< class Value >
	struct coef_b_type : public boost::array< Value , 1 >
	{
		coef_b_type( void )
		{
			(*this)[0] = 1.0;
		}
	};

} // namespace symplectic_euler_coef
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
	class AdjustSizePolicy = adjust_size_initially_tag
	>
class symplectic_euler :
	public symplectic_nystroem_stepper_base
	<
		1 ,
		symplectic_euler< Coor , Momentum , Value , CoorDeriv , MomentumDeriv , Time , Algebra , Operations , AdjustSizePolicy > ,
		Coor , Momentum , Value , CoorDeriv , MomentumDeriv , Time , Algebra , Operations , AdjustSizePolicy
	>
{
public:

		BOOST_ODEINT_SYMPLECTIC_NYSTROEM_STEPPER_TYPEDEFS( symplectic_euler , 1 );

		symplectic_euler( void )
		: stepper_base_type( detail::symplectic_euler_coef::coef_a_type< value_type >() , detail::symplectic_euler_coef::coef_b_type< value_type >() )
		{
		}
};



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_STEPPER_SYMPLECTIC_EULER_HPP_ */
