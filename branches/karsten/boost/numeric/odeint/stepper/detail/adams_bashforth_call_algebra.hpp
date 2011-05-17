/*
 * adams_bashforth_coefficients.hpp
 *
 *  Created on: May 15, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_ADAMS_BASHFORTH_CALL_ALGEBRA_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_ADAMS_BASHFORTH_CALL_ALGEBRA_HPP_

#include <boost/array.hpp>


namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< size_t Step , class Algebra , class Operations >
struct call_algebra;

template< class Algebra , class Operations >
struct call_algebra< 1 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Value >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Value &dt ) const
	{

	}
};


template< class Algebra , class Operations >
struct call_algebra< 2 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Value >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Value &dt ) const
	{

	}
};


template< class Algebra , class Operations >
struct call_algebra< 3 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Value >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Value &dt ) const
	{

	}
};


template< class Algebra , class Operations >
struct call_algebra< 4 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Value >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Value &dt ) const
	{

	}
};


template< class Algebra , class Operations >
struct call_algebra< 5 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Value >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Value &dt ) const
	{

	}
};


template< class Algebra , class Operations >
struct call_algebra< 6 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Value >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Value &dt ) const
	{

	}
};


template< class Algebra , class Operations >
struct call_algebra< 7 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Value >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Value &dt ) const
	{

	}
};


template< class Algebra , class Operations >
struct call_algebra< 8 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Value >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Value &dt ) const
	{

	}
};




} // detail
} // odeint
} // numeric
} // boost



#endif /* BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_ADAMS_BASHFORTH_CALL_ALGEBRA_HPP_ */
