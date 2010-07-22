/*
 boost header: NUMERIC_ODEINT_STEPPER/controlled_error_stepper.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_ERROR_STEPPER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_ERROR_STEPPER_HPP_INCLUDED

#include <boost/numeric/odeint/algebra/standard_algebra.hpp>
#include <boost/numeric/odeint/algebra/standard_operations.hpp>

#include <boost/numeric/odeint/stepper/explicit_stepper_base.hpp>
#include <boost/numeric/odeint/stepper/detail/macros.hpp>

namespace boost {
namespace numeric {
namespace odeint {

typedef enum
{
    success_step_size_unchanged ,
    step_size_decreased ,
    success_step_size_increased
} controlled_step_result;


template<
	class ErrorStepper ,


	>
class controlled_error_stepper
{
public:

	typedef ErrorStepper error_stepper_type;
	typedef typename error_stepper_type::state_type state_type;
	typedef typename error_stepper_type::time_type time_type;
	typedef typename error_stepper_type::order_type order_type;



	template< class System >
	controlled_step_result try_step( System &sys , state_type &x , const state_type &dxdt , time_type &t , time_type &dt )
	{
		using std::max;

		// adjust size

		m_error_checker.fill_scale( x , dxdt , dt , m_x_scale );

		m_x_tmp = x;
		m_stepper.do_step( system , x , dxdt , t , dt , m_x_err );

		time_type max_rel_err = m_error_checker.get_max_error_ratio( m_x_err , m_x_scale );

		if( max_rel_err > 1.1 )
		{
			// error too large - decrease dt
			// limit scaling factor to 0.2
			dt *= std::max( 0.9 * pow( max_rel_err , -1.0/(m_stepper.order_error()-1.0) ),
					0.2 );

			// reset state
			x = m_x_tmp;
			return step_size_decreased;
		}
		else
		{
			if( max_rel_err < 0.5 )
			{
				//error too small - increase dt and keep the evolution
				t += dt;
				// limit scaling factor to 5.0
				dt *= std::min( 0.9*pow(max_rel_err , -1.0/m_stepper.order_error_step()), 5.0 );
				return success_step_size_increased;
			}
			else
			{
				t += dt;
				return success_step_size_unchanged;
			}
		}
	}

	template< class System >
	controlled_step_result try_step( System &sys , state_type &x , time_type &t , time_type &dt )
	{

        system( x , m_dxdt , t );
        return try_step( system , x , m_dxdt , t , dt );
	}


private:

	time_type m_eps_abs;
	time_type m_eps_rel;
	time_type m_a_x;
	time_type m_a_dxdt;

	state_type m_dxdt;
	state_type m_x_tmp;
	state_type m_x_err;
    state_type m_x_scale;

};

//template<
//    class State ,
//    class Time = double ,
//	class Algebra = standard_algebra< State > ,
//	class Operations = standard_operations< Time > ,
//	class AdjustSizePolicy = adjust_size_initially_tag
//	>
//class explicit_euler
//: public explicit_stepper_base<
//	  explicit_euler< State , Time , Algebra , Operations , AdjustSizePolicy > ,
//	  1 , State , Time , Algebra , Operations , AdjustSizePolicy >
//{
//public :
//
//	BOOST_ODEINT_EXPLICIT_STEPPERS_TYPEDEFS( explicit_euler , 1 );
//
//	template< class System >
//	void do_step_impl( System system , state_type &x , const state_type &dxdt , time_type t , time_type dt )
//	{
//		algebra_type::for_each2( x , dxdt , typename operations_type::increment1( dt ) );
//	}
//};




} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_STEPPER_CONTROLLED_ERROR_STEPPER_HPP_INCLUDED
