/*
 * adams_bashforth_coefficients.hpp
 *
 *  Created on: May 15, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_ADAMS_BASHFORTH_CALL_ALGEBRA_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_ADAMS_BASHFORTH_CALL_ALGEBRA_HPP_

#include <cassert>


namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< size_t Step , class Algebra , class Operations >
struct adams_bashforth_call_algebra;

template< class Algebra , class Operations >
struct adams_bashforth_call_algebra< 1 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Time >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Time &dt ) const
	{
		typedef typename Coefficients::value_type value_type;
		Algebra::for_each3( out , in , steps[0] , typename Operations::template scale_sum2< value_type , Time >( 1.0 , dt * coef[0] ) );
	}
};


template< class Algebra , class Operations >
struct adams_bashforth_call_algebra< 2 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Time >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Time &dt ) const
	{
		typedef typename Coefficients::value_type value_type;
		Algebra::for_each4( out , in , steps[0] , steps[1] ,
				typename Operations::template scale_sum3< value_type , Time , Time >( 1.0 , dt * coef[0] , dt * coef[1] ) );
	}
};


template< class Algebra , class Operations >
struct adams_bashforth_call_algebra< 3 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Time >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Time &dt ) const
	{
		typedef typename Coefficients::value_type value_type;
		Algebra::for_each5( out , in , steps[0] , steps[1] , steps[2] ,
				typename Operations::template scale_sum4< value_type , Time , Time >( 1.0 , dt * coef[0] , dt * coef[1] , dt * coef[2] ) );
	}
};


template< class Algebra , class Operations >
struct adams_bashforth_call_algebra< 4 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Time >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Time &dt ) const
	{
		typedef typename Coefficients::value_type value_type;
		Algebra::for_each6( out , in , steps[0] , steps[1] , steps[2] , steps[3] ,
				typename Operations::template scale_sum5< value_type , Time , Time , Time >(
						1.0 , dt * coef[0] , dt * coef[1] , dt * coef[2] , dt * coef[3] ) );
	}
};


template< class Algebra , class Operations >
struct adams_bashforth_call_algebra< 5 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Time >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Time &dt ) const
	{
		typedef typename Coefficients::value_type value_type;
		Algebra::for_each7( out , in , steps[0] , steps[1] , steps[2] , steps[3] , steps[4] ,
				typename Operations::template scale_sum6< value_type , Time , Time , Time , Time >(
						1.0 , dt * coef[0] , dt * coef[1] , dt * coef[2] , dt * coef[3] , dt * coef[4] ) );
	}
};


template< class Algebra , class Operations >
struct adams_bashforth_call_algebra< 6 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Time >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Time &dt ) const
	{
		typedef typename Coefficients::value_type value_type;
		Algebra::for_each8( out , in , steps[0] , steps[1] , steps[2] , steps[3] , steps[4] , steps[5] ,
				typename Operations::template scale_sum7< value_type , Time , Time , Time , Time , Time >(
						1.0 , dt * coef[0] , dt * coef[1] , dt * coef[2] , dt * coef[3] , dt * coef[4] , dt * coef[5] ) );
	}
};


template< class Algebra , class Operations >
struct adams_bashforth_call_algebra< 7 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Time >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Time &dt ) const
	{
		assert( false ); // not implemented
//		typedef typename Coefficients::value_type value_type;
//		Algebra::for_each9( out , in , steps[0] , steps[1] , steps[2] , steps[3] , steps[4] , steps[5] , steps[6]
//				typename Operations::template scale_sum8< value_type , Time , Time , Time , Time , Time , Time >(
//						1.0 , dt * coef[0] , dt * coef[1] , dt * coef[2] , dt * coef[3] , dt * coef[4] , dt * coef[5] , dt * coef[6] ) );
	}
};


template< class Algebra , class Operations >
struct adams_bashforth_call_algebra< 8 , Algebra , Operations >
{
	template< class StateIn , class StateOut , class StepStorage , class Coefficients , class Time >
	void operator()( const StateIn &in , StateOut &out , const StepStorage &steps , const Coefficients &coef , const Time &dt ) const
	{
		assert( false ); // not implemented
//		typedef typename Coefficients::value_type value_type;
//		Algebra::for_each10( out , in , steps[0] , steps[1] , steps[2] , steps[3] , steps[4] , steps[5] , steps[6] , steps[7] ,
//				typename Operations::template scale_sum9< value_type , Time , Time , Time , Time , Time , Time , Time >(
//						1.0 , dt * coef[0] , dt * coef[1] , dt * coef[2] , dt * coef[3] , dt * coef[4] , dt * coef[5] , dt * coef[6] , dt * coef[7] ) );
	}
};




} // detail
} // odeint
} // numeric
} // boost



#endif /* BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_ADAMS_BASHFORTH_CALL_ALGEBRA_HPP_ */
