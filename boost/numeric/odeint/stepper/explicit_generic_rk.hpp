/*
 * explicit_generic_rk.hpp
 *
 *  Created on: May 19th, 2011
 *      Author: mario
 */

#ifndef EXPLICIT_GENERIC_RK_HPP_
#define EXPLICIT_GENERIC_RK_HPP_

#include <boost/mpl/vector.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/size_t.hpp>

#include <boost/fusion/container.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/iterator.hpp>
#include <boost/fusion/mpl.hpp>
#include <boost/fusion/sequence.hpp>

#include <boost/array.hpp>

#include <boost/ref.hpp>

#include <boost/numeric/odeint/stepper/base/explicit_stepper_base.hpp>
#include <boost/numeric/odeint/stepper/detail/macros.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/stepper/detail/generic_rk_call_algebra.hpp>
#include <boost/numeric/odeint/stepper/detail/generic_rk_operations.hpp>
//#include "fusion_foreach_performance.hpp"

#include <iostream>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;


namespace boost {
namespace numeric {
namespace odeint {

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
    class Value = double ,
    class Deriv = State ,
    class Time = Value ,
    class Algebra = range_algebra ,
    class Operations = default_operations ,
    class AdjustSizePolicy = adjust_size_initially_tag
    >
class explicit_generic_rk;

struct stage_vector;

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
Order , State , Value , Deriv , Time , Algebra , Operations , AdjustSizePolicy > stepper_base_type;

	typedef typename stepper_base_type::state_type state_type;
	typedef typename stepper_base_type::value_type value_type;
	typedef typename stepper_base_type::deriv_type deriv_type;
	typedef typename stepper_base_type::time_type time_type;
	typedef typename stepper_base_type::algebra_type algebra_type;
	typedef typename stepper_base_type::operations_type operations_type;
	typedef typename stepper_base_type::adjust_size_policy adjust_size_policy;
	typedef typename stepper_base_type::stepper_type stepper_type;


	typedef mpl::range_c< size_t , 1 , StageCount > stage_indices;

    typedef typename fusion::result_of::as_vector
    <
        typename mpl::copy
        <
            stage_indices ,
            mpl::inserter
            <
                mpl::vector0< > ,
                mpl::push_back< mpl::_1 , array_wrapper< double , mpl::_2 > >
            >
        >::type
    >::type coef_a_type;

    typedef boost::array< double , StageCount > coef_b_type;
    typedef boost::array< double , StageCount > coef_c_type;

    typedef typename fusion::result_of::as_vector
    <
        typename mpl::push_back
        <
            typename mpl::copy
            <
                stage_indices,
                mpl::inserter
                <
                    mpl::vector0<> ,
                    mpl::push_back< mpl::_1 , stage_wrapper< Value , mpl::_2 > >
                >
            >::type ,
            stage< Value , StageCount >
        >::type
    >::type stage_vector_base;


    struct stage_vector : public stage_vector_base
    {
        struct do_insertion
        {
            stage_vector_base &m_base;
            const coef_a_type &m_a;
            const coef_c_type &m_c;

            do_insertion( stage_vector_base &base , const coef_a_type &a , const coef_c_type &c )
            : m_base( base ) , m_a( a ) , m_c( c ) { }

            template< class Index >
            void operator()( Index ) const
            {
                //fusion::at< Index >( m_base ) = stage< double , Index::value+1 , intermediate_stage >( m_c[ Index::value ] , fusion::at< Index >( m_a ) );
                fusion::at< Index >( m_base ).c  = m_c[ Index::value ];
                fusion::at< Index >( m_base ).a = fusion::at< Index >( m_a );
            }
        };


        struct do_insertion_from_stage
        {
            stage_vector_base &m_base;
            const stage_vector_base &m_source;

            do_insertion_from_stage( stage_vector_base &base, const stage_vector_base &source )
                    : m_base(base) , m_source( source )
            { }

            template<class Index>
            void operator()(Index) const {
                //fusion::at< Index >( m_base ) = stage< double , Index::value+1 , intermediate_stage >( m_c[ Index::value ] , fusion::at< Index >( m_a ) );
                fusion::at<Index>(m_base).c = fusion::at<Index>(m_source).c;
                fusion::at<Index>(m_base).a = fusion::at<Index>(m_source).a;
            }
        };

        struct print_butcher
        {
            const stage_vector_base &m_base;
            std::ostream &m_os;

            print_butcher( const stage_vector_base &base , std::ostream &os )
                : m_base( base ) , m_os( os )
            { }

            template<class Index>
            void operator()(Index) const {
                m_os << fusion::at<Index>(m_base).c << " | ";
                for( size_t i=0 ; i<Index::value ; ++i )
                    m_os << fusion::at<Index>(m_base).a[i] << " ";
                m_os << std::endl;
            }
        };


        stage_vector( const coef_a_type &a , const coef_b_type &b , const coef_c_type &c )
        {
            typedef mpl::range_c< size_t , 0 , StageCount-1 > indices;
            mpl::for_each< indices >( do_insertion( *this , a , c ) );
            fusion::at_c< StageCount - 1 >( *this ).c = c[ StageCount - 1 ];
            fusion::at_c< StageCount - 1 >( *this ).a = b;
        }

        stage_vector( const stage_vector &s )
        {
            typedef mpl::range_c< size_t , 0 , StageCount > indices;
            mpl::for_each< indices >( do_insertion_from_stage( *this , s ) );
        }

        void print( std::ostream &os ) const
        {
            typedef mpl::range_c< size_t , 0 , StageCount > indices;
            mpl::for_each< indices >( print_butcher( *this , os ) );
        }
    };



    template< class System , class StateIn , class DerivIn , class StateOut >
    struct calculate_stage
    {
        System &system;
        const StateIn &x;
		state_type &x_tmp;
		StateOut &x_out;
		const DerivIn &dxdt;
        state_type *F;
        const double t;
        const double dt;

        calculate_stage( System &_system , const StateIn &_x , const DerivIn &_dxdt , StateOut &_out , 
			state_type &_x_tmp , state_type *_F , const time_type &_t , const time_type &_dt )
        : system( _system ) , x( _x ) , x_tmp( _x_tmp ) , x_out( _out) , dxdt( _dxdt ) , F( _F ) , t( _t ) , dt( _dt )
        {}


        template< typename T , size_t stage_number >
        void inline operator()( stage< T , stage_number > const &stage ) const
        //typename stage_fusion_wrapper< T , mpl::size_t< stage_number > , intermediate_stage >::type const &stage ) const
        {
            if( stage_number > 1 )
            {
				#ifdef BOOST_MSVC
				#pragma warning( disable : 4307 34 )
				#endif
                system( x_tmp , F[stage_number-2] , t + stage.c * dt );
				#ifdef BOOST_MSVC
				#pragma warning( default : 4307 34 )
				#endif
            }
			//std::cout << stage_number-2 << ", t': " << t + stage.c * dt << std::endl;

			if( stage_number < StageCount )
			    detail::template generic_rk_call_algebra< stage_number , algebra_type >()( x_tmp , x , dxdt , F ,
                            detail::generic_rk_scale_sum< stage_number , operations_type , time_type >( stage.a , dt) );
//				    algebra_type::template for_eachn<stage_number>( x_tmp , x , dxdt , F ,
//				            typename operations_type::template scale_sumn< stage_number , time_type >( stage.a , dt ) );
			else
			    detail::template generic_rk_call_algebra< stage_number , algebra_type >()( x_out , x , dxdt , F ,
                            detail::generic_rk_scale_sum< stage_number , operations_type , time_type >( stage.a , dt) );
//                algebra_type::template for_eachn<stage_number>( x_out , x , dxdt , F ,
//                            typename operations_type::template scale_sumn< stage_number , time_type >( stage.a , dt ) );
        }

    };


private:

    void initialize( void )
    {
        boost::numeric::odeint::construct( m_x_tmp );
        for( size_t i = 0 ; i < StageCount-1 ; ++i )
        {
            boost::numeric::odeint::construct( m_F[i] );
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

    static const size_t stage_count = StageCount;

    explicit_generic_rk( const coef_a_type &a , const coef_b_type &b , const coef_c_type &c )
        : m_stages( a , b , c ) , m_x_tmp()

    {
        initialize();
    }

    explicit_generic_rk( const explicit_generic_rk &rk )
        : stepper_base_type( rk ) , m_stages( rk.m_stages ) , m_x_tmp()
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
	void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt , const time_type &t , StateOut &out , const time_type &dt )
	{
        fusion::for_each( m_stages , calculate_stage< System , StateIn , DerivIn , StateOut >
			( system , in , dxdt , out , m_x_tmp , m_F , t , dt ) );
    }

    friend std::ostream& operator << <>( std::ostream &os , const explicit_generic_rk &rk );

private:

    const stage_vector m_stages;
    state_type m_x_tmp;

protected:
    deriv_type m_F[StageCount-1];

};

}
}
}
#endif /* EXPLICIT_GENERIC_RK_HPP_ */
