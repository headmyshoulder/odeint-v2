/*
 * fusion_stepper.hpp
 *
 *  Created on: Nov 21, 2010
 *      Author: mario
 */

#ifndef FUSION_STEPPER_HPP_
#define FUSION_STEPPER_HPP_

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

#include <boost/fusion/container.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/include/mpl.hpp>

#include <boost/array.hpp>

#include <typeinfo>


#include "algebra.hpp"

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

using namespace std;

template< typename State , size_t stage >
class stepper_stage
{

public:

    typedef double value_type;
    typedef State state_type;
    typedef vector< vector< double > > parameter_array_type2D;
    typedef vector< double > parameter_array_type1D;


    void init( const parameter_array_type2D &a , const parameter_array_type1D &c )
    {
        copy( a[stage-1].begin() , a[stage-1].end() , &m_a_row[0] );
        m_c = c[stage-1];
    }

    template< typename System >
    void operator() ( System &system , state_type &x_tmp , const state_type &x , state_type *k_vector ,
    					double t , const double dt ) const
    {
        system( x_tmp , k_vector[stage-1] , t + m_c * dt );
        algebra<state_type , stage>::foreach< stage >( x_tmp , x , m_a_row , k_vector , dt);
    }

private:

    boost::array< value_type , stage > m_a_row;
    value_type m_c;

};

template< typename State >
class stepper_stage< State , 1 >
{

public:

    typedef double value_type;
    typedef State state_type;
    typedef vector< vector< double > > parameter_array_type2D;
    typedef vector< double > parameter_array_type1D;


    void init( const parameter_array_type2D &a , const parameter_array_type1D &c )
    {
    	copy( a[0].begin() , a[0].end() , &m_a_row[0] );
        m_c = c[0];
    }

    template< typename System >
    void operator() ( System &system , state_type &x_tmp , const state_type &x , state_type *k_vector ,
    					double t , const double dt ) const
    {
        system( x , k_vector[0] , t + m_c * dt );
        algebra<state_type , 1>::foreach( x_tmp , x , m_a_row , k_vector , dt);
    }

private:

    boost::array< value_type , 1 > m_a_row;
    value_type m_c;

};


template< typename State , size_t stage >
class stepper_last_stage
{

public:

    typedef double value_type;
    typedef State state_type;
    typedef vector< vector< double > > parameter_array_type2D;
    typedef vector< double > parameter_array_type1D;


    void init( const parameter_array_type1D &b , const parameter_array_type1D &c )
    {
        copy( b.begin() , b.end() , &m_b[0] );
        m_c = c[stage-1];
    }

    template< typename System >
    void operator() ( System &system , state_type &x_tmp , state_type &x , state_type *k_vector ,
                        double t , const double dt ) const
    {
        system( x_tmp , k_vector[stage-1] , t + m_c * dt );
        algebra<state_type , stage>::foreach( x , x , m_b , k_vector , dt);
    }

private:

    boost::array< value_type , stage > m_b;
    value_type m_c;

};


template< typename State >
class stepper_last_stage< State , 1 >
{

public:

    typedef double value_type;
    typedef State state_type;
    typedef vector< vector< double > > parameter_array_type2D;
    typedef vector< double > parameter_array_type1D;


    void init( const parameter_array_type1D &b , const parameter_array_type1D &c )
    {
        copy( b.begin() , b.end() , &m_b[0] );
        m_c = c[0];
    }

    template< typename System >
    void operator() ( System &system , state_type &x_tmp , state_type &x , state_type k_vector[1] ,
                        double t , const double dt ) const
    {
        system( x , k_vector[0] , t + m_c * dt );
        algebra<state_type , 1>::foreach( x , x , m_b , k_vector , dt);
    }

private:

    boost::array< value_type , 1 > m_b;
    value_type m_c;

};








template< class T , class Constant >
struct stepper_stage_wrapper
{
	typedef stepper_stage< T , Constant::value > type;
};




/* Runge Kutta Stepper - consisting of several stages */

template< typename State , size_t stage_count >
class runge_kutta_stepper
{

	typedef State state_type;


	typedef typename fusion::result_of::as_vector
    <
        typename mpl::copy
        <
            mpl::range_c< size_t , 1 , N > ,
            mpl::inserter
            <
                mpl::vector0<> ,
                mpl::push_back< mpl::_1 , array_wrapper< double , mpl::_2 > >
            >
        >::type
    >::type coef_a_type;
	typedef boost::array< double , N > coef_b_type;
	typedef boost::array< double , N > coef_c_type;


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

    	template< class Stage >
    	void operator()( Stage const &stage ) const
    	{
    		stage( system , x_tmp , x , k_vector , t , dt );
    	}

    };

    struct init_stage {

    	const parameter_array_type2D &m_a;
    	const parameter_array_type1D &m_c;

    	init_stage( const parameter_array_type2D &a , const parameter_array_type1D &c )
    		: m_a( a ) , m_c( c )
    	{ }

    	template< typename Stage >
    	void operator() ( Stage &stage ) const
    	{
    		stage.init( m_a , m_c );
    	}
    };

public:

    typedef mpl::range_c< size_t , 1 , stage_count > indices;
    typedef typename fusion::result_of::as_vector <
    	typename mpl::copy // create mpl::vector< stepper_stage< 0 >, stepper_stage< 1 > , .. stepper_stage< stage_count-1 > >
    	<
    		indices ,
    		mpl::inserter
    		<
    			mpl::vector0<> ,
    			mpl::push_back< mpl::_1 , stepper_stage_wrapper< State , mpl::_2 > >
      	  	>
    	>::type
    >::type stage_vector;
    typedef stepper_last_stage< State , stage_count > last_stage_type;

    runge_kutta_stepper( const coef_a_type &a ,
							const coef_b_type &b ,
							const coef_c_type &c )
    {
    	fusion::for_each( m_stages , init_stage( a , c ) );
    	m_last_stage.init( b , c );
    }

    template< class System >
    void do_step( System &system , state_type &x , double t , const double dt )
    {
        fusion::for_each( m_stages , calculate_stage< System >( system , x , m_x_tmp , m_k_vector , t , dt ) );
        m_last_stage( system , m_x_tmp , x , m_k_vector , t , dt );
    }


private:

	state_type m_k_vector[stage_count];
	state_type m_x_tmp;
	stage_vector m_stages;
	last_stage_type m_last_stage;
};

#endif /* FUSION_STEPPER_HPP_ */
