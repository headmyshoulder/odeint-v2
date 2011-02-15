/*
 * symplectic_nystroem.hpp
 *
 *  Created on: Feb 12, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_SYMPLECTIC_NYSTROEM_STEPPER_BASE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_SYMPLECTIC_NYSTROEM_STEPPER_BASE_HPP_INCLUDED

#include <boost/ref.hpp>
#include <boost/array.hpp>

#include <boost/numeric/odeint/util/size_adjuster.hpp>
#include <boost/numeric/odeint/util/construct.hpp>
#include <boost/numeric/odeint/util/destruct.hpp>
#include <boost/numeric/odeint/util/copy.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

#include <boost/numeric/odeint/stepper/detail/is_pair.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template<
	size_t NumOfStages ,
	class Stepper ,
	class Coor ,
	class Momentum ,
	class Value ,
	class CoorDeriv ,
	class MomentumDeriv ,
	class Time ,
	class Algebra ,
	class Operations ,
	class AdjustSizePolicy
	>
class symplectic_nystroem_stepper_base
{

	void initialize( void )
	{
		boost::numeric::odeint::construct( m_dqdt );
		boost::numeric::odeint::construct( m_dpdt );
		m_coor_deriv_adjuster.register_state( 0 , m_dqdt );
		m_momentum_deriv_adjuster.register_state( 0 , m_dpdt );
	}

	void copy( const symplectic_nystroem_stepper_base &b )
	{
		boost::numeric::odeint::copy( b.m_dqdt , m_dqdt );
		boost::numeric::odeint::copy( b.m_dpdt , m_dpdt );
	}

public:

	const static size_t num_of_stages = NumOfStages;
	typedef Coor coor_type;
	typedef Momentum momentum_type;
	typedef std::pair< coor_type , momentum_type > state_type;
	typedef CoorDeriv coor_deriv_type;
	typedef MomentumDeriv momentum_deriv_type;
	typedef std::pair< coor_deriv_type , momentum_deriv_type > deriv_type;
	typedef Value value_type;
	typedef Time time_type;
	typedef Algebra algebra_type;
	typedef Operations operations_type;
	typedef AdjustSizePolicy adjust_size_policy;
	typedef Stepper stepper_type;
	typedef stepper_tag stepper_category;

	typedef boost::array< value_type , num_of_stages > coef_type;

	symplectic_nystroem_stepper_base( const coef_type &coef_a , const coef_type &coef_b )
	: m_coef_a( coef_a ) , m_coef_b( coef_b ) ,
	  m_coor_deriv_adjuster() , m_momentum_deriv_adjuster() ,
	  m_dqdt() , m_dpdt()
	{
		initialize();
	}

	symplectic_nystroem_stepper_base( const symplectic_nystroem_stepper_base &b )
	: m_coef_a( b.m_coef_a ) , m_coef_b( b.m_coef_b ) ,
	  m_coor_deriv_adjuster() , m_momentum_deriv_adjuster(),
	  m_dqdt() , m_dpdt()
	{
		initialize();
		copy( b );
	}

	~symplectic_nystroem_stepper_base( void )
	{
		boost::numeric::odeint::destruct( m_dqdt );
		boost::numeric::odeint::destruct( m_dpdt );
	}

	symplectic_nystroem_stepper_base& operator=( const symplectic_nystroem_stepper_base &b )
	{
		copy( b );
		return *this;
	}


	template< class System , class StateIn , class StateOut >
	void do_step( System system , const StateIn &in , const time_type &t , const StateOut &out , const time_type &dt )
	{
		do_step( system , in , t , out , dt );
	}


	template< class System , class StateIn , class StateOut >
	void do_step( System system , const StateIn &in , const time_type &t , StateOut &out , const time_type &dt )
	{
		typedef typename boost::unwrap_reference< System >::type system_type;
		do_step_impl( system , in , t , out , dt , typename detail::is_pair< system_type >::type() );
	}




	template< class System , class StateInOut >
	void do_step( System system , const StateInOut &state , const time_type &t , const time_type &dt )
	{
		do_step( system , state , t , dt );
	}

	template< class System , class StateInOut >
	void do_step( System system , StateInOut &state , const time_type &t , const time_type &dt )
	{
		typedef typename boost::unwrap_reference< System >::type system_type;
		do_step_impl( system , state , t , state , dt , typename detail::is_pair< system_type >::type() );
	}

	template< class StateType >
	void adjust_size( const StateType &x )
	{
		m_coor_deriv_adjuster.adjust_size( x );
		m_momentum_deriv_adjuster.adjust_size( x );
	}

	const coef_type& coef_a( void ) const { return m_coef_a; }
	const coef_type& coef_b( void ) const { return m_coef_b; }

private:

	// stepper for systems with function for dq/dt = f(p) and dp/dt = -f(q)
	template< class System , class StateIn , class StateOut >
	void do_step_impl( System system , const StateIn &in , const time_type &t , StateOut &out , const time_type &dt , boost::mpl::true_ )
	{
    	typedef typename boost::unwrap_reference< System >::type system_type;
    	typedef typename boost::unwrap_reference< typename system_type::first_type >::type coor_deriv_func_type;
    	typedef typename boost::unwrap_reference< typename system_type::second_type >::type momentum_deriv_func_type;
    	system_type &sys = system;
    	coor_deriv_func_type &coor_func = sys.first;
    	momentum_deriv_func_type &momentum_func = sys.second;


		m_coor_deriv_adjuster.adjust_size_by_policy( in , adjust_size_policy() );
		m_momentum_deriv_adjuster.adjust_size_by_policy( in , adjust_size_policy() );

		// ToDo: check sizes?

		for( size_t l=0 ; l<num_of_stages ; ++l )
		{
			if( l == 0 )
			{
				coor_func( in.second , m_dqdt );
				typename algebra_type::for_each2()( out.first , in.first , m_dqdt ,
					typename operations_type::template scale_sum2< value_type , time_type >( 1.0 , coef_a[l] * dt ) );
				momentum_func( out.first , m_dqdt );
				typename algebra_type::for_each2()( out.second , in.second , m_dqdt ,
					typename operations_type::template scale_sum2< value_type , time_type >( 1.0 , coef_b[l] * dt ) );
			}
			else
			{
				coor_func( out.second , m_dqdt );
				typename algebra_type::for_each2()( out.first , out.first , m_dqdt ,
					typename operations_type::template scale_sum2< value_type , time_type >( 1.0 , coef_a[l] * dt ) );
				momentum_func( out.first , m_dqdt );
				typename algebra_type::for_each2()( out.second , out.second , m_dqdt ,
					typename operations_type::template scale_sum2< value_type , time_type >( 1.0 , coef_b[l] * dt ) );
			}
		}
	}


	// stepper for systems with only function dp /dt = -f(q), dq/dt = p
	template< class System , class StateIn , class StateOut >
	void do_step_impl( System system , const StateIn &in , const time_type &t , StateOut &out , const time_type &dt , boost::mpl::false_ )
	{
    	typedef typename boost::unwrap_reference< System >::type momentum_deriv_func_type;
    	momentum_deriv_func_type &momentum_func = system;


		m_coor_deriv_adjuster.adjust_size_by_policy( in , adjust_size_policy() );
		m_momentum_deriv_adjuster.adjust_size_by_policy( in , adjust_size_policy() );

		// ToDo: check sizes?


		for( size_t l=0 ; l<num_of_stages ; ++l )
		{
			if( l == 0 )
			{
				typename algebra_type::for_each3()( out.first , in.first , in.second ,
					typename operations_type::template scale_sum2< value_type , time_type >( 1.0 , m_coef_a[l] * dt ) );
				momentum_func( out.first , m_dqdt );
				typename algebra_type::for_each3()( out.second , in.second , m_dqdt ,
					typename operations_type::template scale_sum2< value_type , time_type >( 1.0 , m_coef_b[l] * dt ) );
			}
			else
			{
				typename algebra_type::for_each3()( out.first , out.first , out.second ,
					typename operations_type::template scale_sum2< value_type , time_type >( 1.0 , m_coef_a[l] * dt ) );
				momentum_func( out.first , m_dqdt );
				typename algebra_type::for_each3()( out.second , out.second , m_dqdt ,
					typename operations_type::template scale_sum2< value_type , time_type >( 1.0 , m_coef_b[l] * dt ) );
			}
		}
	}

	// ToDo : check if a and b are appropriate names for the coefficients
	const coef_type m_coef_a;
	const coef_type m_coef_b;

	size_adjuster< coor_deriv_type , 1 > m_coor_deriv_adjuster;
	size_adjuster< momentum_deriv_type , 1 > m_momentum_deriv_adjuster;
	coor_deriv_type m_dqdt;
	momentum_deriv_type m_dpdt;
};

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_SYMPLECTIC_NYSTROEM_STEPPER_BASE_HPP_INCLUDED
