#include <boost/mpl/vector.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/size_t.hpp>

#include <boost/fusion/container.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/adapted.hpp>


#include <algorithm>
#include <iostream>
#include <string>

#include <boost/array.hpp>
#include <typeinfo>


#include "fusion_algebra.hpp"

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

using namespace std;


struct first_stage {};
struct intermediate_stage {};
struct last_stage {};

ostream& operator<<( ostream& out , const first_stage& )
{
	out << "first stage";
	return out;
}

ostream& operator<<( ostream& out , const intermediate_stage& )
{
	out << "intermediate stage";
	return out;
}

ostream& operator<<( ostream& out , const last_stage& )
{
	out << "last stage";
	return out;
}



template< class T , class Constant >
struct array_wrapper
{
    typedef typename boost::array< T , Constant::value > type;
};

template< class T , class Constant , class StageCategory >
struct stage_fusion_wrapper
{
    typedef typename fusion::vector< size_t , T , boost::array< T , Constant::value > , StageCategory > type;
};

struct print_stage
{
    template< typename Stage >
    void operator() ( const Stage &stage ) const
    {
        std::cout<< "stage: " << fusion::at_c<0>( stage ) << " c:" << fusion::at_c<1>( stage ) << " ";
        for( size_t i=0 ; i < fusion::at_c<2>( stage ).size() ; ++i )
            std::cout << "a["<<i<<"]:" <<  fusion::at_c<2>( stage )[i] << " ";
        std::cout << std::endl;
    }
};

struct print_sequence
{
	size_t m_indent;

	print_sequence( size_t indent = 0 )
	: m_indent( indent ) { }

	template< typename T >
	void operator()( T ) const
	{
		string in;
		for( size_t i=0 ; i<m_indent ; ++i ) in += "  ";
		cout << in << typeid(T).name() << endl;
	}
};

struct print_values
{
	size_t m_indent;

	print_values( size_t indent = 0 )
	: m_indent( indent ) { }

	template< typename T >
	void operator()( const T &t ) const
	{
		typedef typename fusion::traits::is_sequence< T >::type seq_type;
		print_val( t , seq_type() );
	}

	template< class R >
	void print_val( const R &r , mpl::false_ ) const
	{
		string in;
		for( size_t i=0 ; i<m_indent ; ++i ) in += "  ";
		cout << in << r << endl;
	}

	template< class S >
	void print_val( const S &s , mpl::true_ ) const
	{
		fusion::for_each( s , print_values( m_indent + 1 ) );
	}
};


template< class StateType , size_t stage_count >
class runge_kutta_stepper
{

public:

    typedef StateType state_type;

    typedef mpl::range_c< size_t , 1 , stage_count > stage_indices;

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

    typedef boost::array< double , stage_count > coef_b_type;
    typedef boost::array< double , stage_count > coef_c_type;

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
                    mpl::push_back< mpl::_1 , stage_fusion_wrapper< double , mpl::_2 , intermediate_stage > >
                >
            >::type ,
            typename stage_fusion_wrapper< double , mpl::size_t< stage_count > , last_stage >::type
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
    			fusion::at_c< 0 >( fusion::at< Index >( m_base ) ) = Index::value;
    			fusion::at_c< 1 >( fusion::at< Index >( m_base ) ) = m_c[ Index::value ];
    			fusion::at_c< 2 >( fusion::at< Index >( m_base ) ) = fusion::at< Index >( m_a );
    		}
    	};

    	stage_vector( const coef_a_type &a , const coef_b_type &b , const coef_c_type &c )
    	{
    		typedef mpl::range_c< size_t , 0 , stage_count - 1 > indices;
    		mpl::for_each< indices >( do_insertion( *this , a , c ) );
			fusion::at_c< 0 >( fusion::at_c< stage_count - 1 >( *this ) ) = stage_count - 1 ;
			fusion::at_c< 1 >( fusion::at_c< stage_count - 1 >( *this ) ) = c[ stage_count - 1 ];
			fusion::at_c< 2 >( fusion::at_c< stage_count - 1 >( *this ) ) = b;
    	}
    };



    template< class System >
    struct calculate_stage
    {
        System &system;
        state_type &x , &x_tmp;
        state_type *k_vector;
        const double t;
        const double dt;

        calculate_stage( System &_system , state_type &_x , state_type &_x_tmp , state_type *_k_vector ,
                            const double _t , const double _dt )
        : system( _system ) , x( _x ) , x_tmp( _x_tmp ) , k_vector( _k_vector ) , t( _t ) , dt( _dt )
        {}


        template< typename T , size_t stage_number >
        void operator()( fusion::vector< size_t , T , boost::array< T , stage_number > , intermediate_stage > const &stage ) const
        //typename stage_fusion_wrapper< T , mpl::size_t< stage_number > , intermediate_stage >::type const &stage ) const
        {
            double c = fusion::at_c< 1 >( stage );

            if( stage_number == 1 )
                system( x , k_vector[stage_number-1] , t + c * dt );
            else
                system( x_tmp , k_vector[stage_number-1] , t + c * dt );

            fusion_algebra<stage_number>::foreach( x_tmp , x , fusion::at_c< 2 >( stage ) , k_vector , dt);
        }


        template< typename T , size_t stage_number >
        void operator()( fusion::vector< size_t , T , boost::array< T , stage_number > , last_stage > const &stage ) const
        //void operator()( typename stage_fusion_wrapper< T , mpl::size_t< stage_number > , last_stage >::type const &stage ) const
        {
            double c = fusion::at_c< 1 >( stage );

            if( stage_number == 1 )
                system( x , k_vector[stage_number-1] , t + c * dt );
            else
                system( x_tmp , k_vector[stage_number-1] , t + c * dt );

            fusion_algebra<stage_number>::foreach( x , x , fusion::at_c< 2 >( stage ) , k_vector , dt);
        }


    };




public:

    runge_kutta_stepper( const coef_a_type &a ,
                            const coef_b_type &b ,
                            const coef_c_type &c )
    : m_a( a ) , m_b( b ) , m_c( c ) ,
      m_stages( a , b , c )

    { }


    template< class System >
    void do_step( System &system , state_type &x , double t , const double dt )
    {
        fusion::for_each( m_stages , calculate_stage< System >( system , x , m_x_tmp , m_k_vector , t , dt ) );
    }


    void print_vals()
    {
    	cout << "Typeof state_vector_base : " << endl;
    	mpl::for_each< stage_vector_base >( print_sequence( 1 ) );

    	cout << "m_a : " << endl;
    	fusion::for_each( m_a , print_values( 1 ) );

    	cout << "m_b : " << endl;
    	fusion::for_each( m_b , print_values( 1 ) );

    	cout << "m_c : " << endl;
    	fusion::for_each( m_c , print_values( 1 ) );

    	cout << "m_stages : " << endl;
    	fusion::for_each( m_stages , print_values( 1 ) );
    }



private:

    const coef_a_type m_a;
    const coef_b_type m_b;
    const coef_c_type m_c;
    const stage_vector m_stages;
    state_type m_x_tmp;
    state_type m_k_vector[stage_count];
};
