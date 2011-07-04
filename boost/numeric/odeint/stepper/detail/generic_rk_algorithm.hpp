/*
 * generic_rk_algorithm.hpp
 *
 *  Created on: Jun 3, 2011
 *      Author: mario
 */

#ifndef GENERIC_RK_ALGORITHM_HPP_
#define GENERIC_RK_ALGORITHM_HPP_

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

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/stepper/detail/generic_rk_call_algebra.hpp>
#include <boost/numeric/odeint/stepper/detail/generic_rk_operations.hpp>

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

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
    class Value = double ,
    class Algebra = range_algebra,
    class Operations = default_operations
    >
class generic_rk_algorithm {

public:
    typedef mpl::range_c< size_t , 1 , StageCount > stage_indices;

    typedef typename fusion::result_of::as_vector
    <
        typename mpl::copy
        <
            stage_indices ,
            mpl::inserter
            <
                mpl::vector0< > ,
                mpl::push_back< mpl::_1 , array_wrapper< Value , mpl::_2 > >
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



    template< class System , class StateIn , class StateTemp , class DerivIn , class Deriv ,
               class StateOut , class Time >
    struct calculate_stage
    {
        System &system;
        const StateIn &x;
        StateTemp &x_tmp;
        StateOut &x_out;
        const DerivIn &dxdt;
        Deriv *F;
        const Time t;
        const Time dt;

        calculate_stage( System &_system , const StateIn &_x , const DerivIn &_dxdt , StateOut &_out ,
            StateTemp &_x_tmp , Deriv *_F , const Time &_t , const Time &_dt )
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
                detail::template generic_rk_call_algebra< stage_number , Algebra >()( x_tmp , x , dxdt , F ,
                            detail::generic_rk_scale_sum< stage_number , Operations , Time >( stage.a , dt) );
//                  algebra_type::template for_eachn<stage_number>( x_tmp , x , dxdt , F ,
//                          typename operations_type::template scale_sumn< stage_number , time_type >( stage.a , dt ) );
            else
                detail::template generic_rk_call_algebra< stage_number , Algebra >()( x_out , x , dxdt , F ,
                            detail::generic_rk_scale_sum< stage_number , Operations , Time >( stage.a , dt) );
//                algebra_type::template for_eachn<stage_number>( x_out , x , dxdt , F ,
//                            typename operations_type::template scale_sumn< stage_number , time_type >( stage.a , dt ) );
        }

    };

    generic_rk_algorithm( const coef_a_type &a , const coef_b_type &b , const coef_c_type &c )
        : m_stages( a , b , c )
    { }

    template< class System , class StateIn , class DerivIn , class Time , class StateOut , class StateTemp , class Deriv >
    void inline do_step( System system , const StateIn &in , const DerivIn &dxdt ,
                    const Time &t , StateOut &out , const Time &dt ,
                    StateTemp &x_tmp , Deriv F[StageCount-1] )
    {
        fusion::for_each( m_stages , calculate_stage<
                System , StateIn , StateTemp , DerivIn , Deriv , StateOut , Time >
            ( system , in , dxdt , out , x_tmp , F , t , dt ) );
    }

private:
    stage_vector m_stages;
};


}
}
}
}

#endif /* GENERIC_RK_ALGORITHM_HPP_ */
