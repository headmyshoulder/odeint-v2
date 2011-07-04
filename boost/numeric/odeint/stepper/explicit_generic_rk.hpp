/*
 * explicit_generic_rk.hpp
 *
 *  Created on: May 19th, 2011
 *      Author: mario
 */

#ifndef EXPLICIT_GENERIC_RK_HPP_
#define EXPLICIT_GENERIC_RK_HPP_

#include <boost/array.hpp>

#include <boost/ref.hpp>

#include <boost/numeric/odeint/stepper/base/explicit_stepper_base.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/stepper/detail/generic_rk_algorithm.hpp>

//#include "fusion_foreach_performance.hpp"

#include <iostream>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace boost {
namespace numeric {
namespace odeint {

//forward declarations

template<
    size_t StageCount,
    size_t Order,
    class State ,
    class Value = double ,
    class Deriv = State ,
    class Time = Value ,
    class Algebra = range_algebra ,
    class Operations = default_operations ,
    class AdjustSizePolicy = adjust_size_initially_tag
    >
class explicit_generic_rk;

struct stage_vector;

template< class T , class Constant >
struct array_wrapper
{
    typedef const typename boost::array< T , Constant::value > type;
};

template< class T , size_t i >
struct stage
{
    T c;
    boost::array< T , i > a;
};


template< class T , class Constant >
struct stage_wrapper
{
    typedef stage< T , Constant::value > type;
};


template<
    size_t StageCount,
    size_t Order,
    class State ,
    class Value ,
    class Deriv ,
    class Time ,
    class Algebra ,
    class Operations ,
    class AdjustSizePolicy
    >
std::ostream& operator <<( std::ostream &os ,
        const explicit_generic_rk< StageCount , Order , State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > &rk )
{
      os << "Generic RK with " << rk.stage_count << " stages." << std::endl;
      os << "Butcher Tableau: " << std::endl;
      rk.m_stages.print( os );
      return os;
}


template<
	size_t StageCount,
	size_t Order,
    class State ,
    class Value ,
    class Deriv ,
    class Time ,
	class Algebra ,
	class Operations ,
	class AdjustSizePolicy
	>
class explicit_generic_rk : public explicit_stepper_base<
	  explicit_generic_rk< StageCount , Order , State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > ,
	  Order , State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy >
{

public:

	typedef explicit_stepper_base<
	        explicit_generic_rk< StageCount , Order , State , Value , Deriv ,Time , Algebra , Operations , AdjustSizePolicy > ,
	        Order , State , Value , Deriv , Time , Algebra ,
	        Operations , AdjustSizePolicy > stepper_base_type;

	typedef typename stepper_base_type::state_type state_type;
	typedef typename stepper_base_type::value_type value_type;
	typedef typename stepper_base_type::deriv_type deriv_type;
	typedef typename stepper_base_type::time_type time_type;
	typedef typename stepper_base_type::algebra_type algebra_type;
	typedef typename stepper_base_type::operations_type operations_type;
	typedef typename stepper_base_type::adjust_size_policy adjust_size_policy;
	typedef typename stepper_base_type::stepper_type stepper_type;

	typedef detail::generic_rk_algorithm< StageCount , Value , Algebra , Operations > rk_algorithm_type;

	typedef typename rk_algorithm_type::coef_a_type coef_a_type;
    typedef typename rk_algorithm_type::coef_b_type coef_b_type;
    typedef typename rk_algorithm_type::coef_c_type coef_c_type;

    static const size_t stage_count = StageCount;


private:

    void initialize( void )
    {
        boost::numeric::odeint::construct( m_x_tmp );
        m_state_adjuster.register_state( 0 , m_x_tmp );
        for( size_t i = 0 ; i < StageCount-1 ; ++i )
        {
            boost::numeric::odeint::construct( m_F[i] );
            m_deriv_adjuster.register_state( i , m_F[i] );
        }
    }

    void copy( const explicit_generic_rk &rk )
    {
        boost::numeric::odeint::copy( rk.m_x_tmp , m_x_tmp );
        for( size_t i = 0 ; i < StageCount-1 ; ++i )
        {
            boost::numeric::odeint::copy( rk.m_F[i] , m_F[i] );
        }
    }


public:

    explicit_generic_rk( const coef_a_type &a , const coef_b_type &b , const coef_c_type &c )
        : m_rk_algorithm( a , b , c ) , m_x_tmp()

    {
        initialize();
    }

    explicit_generic_rk( const explicit_generic_rk &rk )
        : stepper_base_type( rk ) , m_rk_algorithm( rk.m_rk_algorithm) , m_x_tmp()
    {
        initialize();
        copy( rk );
    }

    explicit_generic_rk& operator=( const explicit_generic_rk &rk )
    {
        stepper_base_type::operator=( rk );
        copy( rk );
        return *this;
    }

    ~explicit_generic_rk( void )
    {
        boost::numeric::odeint::destruct( m_x_tmp );
        for( size_t i = 0 ; i < StageCount-1 ; ++i )
        {
            boost::numeric::odeint::destruct( m_F[i] );
        }
    }


    template< class System , class StateIn , class DerivIn , class StateOut >
	void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt ,
	                     const time_type &t , StateOut &out , const time_type &dt )
	{
        typedef typename boost::unwrap_reference< System >::type unwrapped_system_type;
        unwrapped_system_type &sys = system;

        m_deriv_adjuster.adjust_size_by_policy( in , adjust_size_policy() );
        m_state_adjuster.adjust_size_by_policy( in , adjust_size_policy() );

        // actual calculation done in generic_rk.hpp
        m_rk_algorithm.do_step( sys , in , dxdt , t , out , dt , m_x_tmp , m_F );
    }

    template< class StateType >
    void adjust_size( const StateType &x )
    {
        m_deriv_adjuster.adjust_size( x );
        m_state_adjuster.adjust_size( x );
        stepper_base_type::adjust_size( x );
    }

    friend std::ostream& operator << <>( std::ostream &os , const explicit_generic_rk &rk );

private:

    rk_algorithm_type m_rk_algorithm;

    size_adjuster< deriv_type , StageCount-1 > m_deriv_adjuster;
    size_adjuster< state_type , 1 > m_state_adjuster;

    state_type m_x_tmp;
    deriv_type m_F[StageCount-1];

};

}
}
}
#endif /* EXPLICIT_GENERIC_RK_HPP_ */
