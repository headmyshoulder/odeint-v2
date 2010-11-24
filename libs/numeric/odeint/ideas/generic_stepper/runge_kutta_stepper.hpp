
#include <boost/mpl/vector.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/zip_view.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/insert_range.hpp>

#include <vector>
#include <algorithm>
#include <iostream>

#include <boost/fusion/container.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/include/mpl.hpp>

#include <boost/array.hpp>

#include <typeinfo>


#include "algebra.hpp"

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

using namespace std;


template< class T , class Constant >
struct array_wrapper
{
    typedef typename boost::array< T , Constant::value > type;
};


template< class T , class Constant >
struct stage_wrapper
{
    typedef typename  fusion::result_of::as_vector< mpl::vector< size_t , T , boost::array< T , Constant::value > > >::type type;
};

struct print_stage
{
    template< typename Stage >
    void operator() ( const Stage &stage ) const
    {
        std::cout<< fusion::at_c<0>( stage ) << " " << fusion::at_c<1>( stage ) << " " << std::endl;
    }
};


/* Runge Kutta Stepper - consisting of several stages */

template< typename State , size_t stage_count >
class runge_kutta_stepper
{
public:
    typedef State state_type;

    typedef mpl::range_c< size_t , 1 , stage_count > stage_indices;

    typedef typename fusion::result_of::as_vector
    <
        typename mpl::copy
        <
            stage_indices ,
            mpl::inserter
            <
                mpl::vector0<> ,
                mpl::push_back< mpl::_1 , array_wrapper< double , mpl::_2 > >
            >
        >::type
    >::type coef_a_type;

    typedef boost::array< double , stage_count > coef_b_type;

    typedef boost::array< double , stage_count > coef_c_type;


    typedef typename fusion::result_of::as_vector<
        typename mpl::copy<
            stage_indices,
            mpl::inserter
              <
                mpl::vector0<> ,
                mpl::push_back< mpl::_1 , stage_wrapper< double , mpl::_2 > >
              >
            >::type
        >::type stage_vector_base;


    struct stage_vector_type : public stage_vector_base
    {
        struct init_stage{

            const coef_a_type &m_a;
            const coef_b_type &m_c;

            init_stage( const coef_a_type &a , const coef_b_type &b , const coef_c_type &c )
                : m_a( a ) , m_c( c ) {}


            template< typename Stage >
            void operator() ( const Stage &stage ) const
            {
                const size_t n = fusion::at_c<2>( stage ).size();
                fusion::at_c<0>( stage ) = n;
                fusion::at_c<1>( stage ) = m_c[n];
                fusion::at_c<2>( stage ) = m_a[n];
            }

        };


        stage_vector_type( const coef_a_type &a , const coef_b_type &b , const coef_c_type &c )
        {
            fusion::for_each( static_cast< stage_vector_base& >( *this ) , init_stage( a , b , c ) );
        }
    };


//    template< class System >
//    struct calculate_stage
//    {
//        System &system;
//        state_type &x , &x_tmp;
//        state_type *k_vector;
//        const double t;
//        const double dt;
//
//        calculate_stage( System &_system , state_type &_x , state_type &_x_tmp , state_type *_k_vector ,
//                            const double _t , const double _dt )
//        : system( _system ) , x( _x ) , x_tmp( _x_tmp ) , k_vector( _k_vector ) , t( _t ) , dt( _dt )
//        {}
//
//        template< class Stage >
//        void operator()( Stage const &stage ) const
//        {
//            stage( system , x_tmp , x , k_vector , t , dt );
//        }
//
//    };


public:

    runge_kutta_stepper( const coef_a_type &a ,
                            const coef_b_type &b ,
                            const coef_c_type &c )
        : m_stages( a , b , c )
    { }

    template< class System >
    void do_step( System &system , state_type &x , double t , const double dt )
    {
//        fusion::for_each( m_stages , calculate_stage< System >( system , x , m_x_tmp , m_k_vector , t , dt ) );
//        m_last_stage( system , m_x_tmp , x , m_k_vector , t , dt );
    }

    void print_vals()
    {
        fusion::for_each( m_stages , print_stage() );
    }


private:

    state_type m_k_vector[stage_count];
    state_type m_x_tmp;
    const stage_vector_type m_stages;
//    last_stage_type m_last_stage;
};
