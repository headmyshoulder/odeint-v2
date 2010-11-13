/*
 * generic_stepper.hpp
 *
 *  Created on: Nov 12, 2010
 *      Author: mario
 */

#ifndef GENERIC_STEPPER_HPP_
#define GENERIC_STEPPER_HPP_


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

#include "algebra.hpp"

namespace mpl = boost::mpl;

using namespace std;


/*
 * c_1 |
 * c_2 | a_{2,1}
 * c_3 | a_{3,1} a_{3,2}
 * ...
 * c_s | a_{s,1} a_{s,2} ... a_{s,s-1}
 * -----------------------------------
 *     | b_1     b_2         b_{s-1}    b_s
 */


/* Runge Kutta Stages */


template< typename State , size_t stage >
class stepper_stage
{

public:

    typedef double value_type;
    typedef State state_type;
    typedef vector< vector< double > > parameter_array_type2D;
    typedef vector< double > parameter_array_type1D;


    void init( parameter_array_type2D &a , parameter_array_type1D &c )
    {
        m_a_row = a[stage-1];
        m_c = c[stage-1];
    }

    template< typename System >
    void operator() ( System &system , state_type &x_tmp , const state_type &x , const state_type *k_vector , double t , const double dt )
    {
        system( x_tmp , k_vector[stage-1] , t + m_c * dt );
        algebra<state_type , stage>::foreach< stage >( x_tmp , x , m_a_row , k_vector , dt);
    }

private:

    vector< value_type > m_a_row;
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


    void init( parameter_array_type2D &a , parameter_array_type1D &c )
    {
    	m_a_row = a[0];
        m_c = c[0];
    }

    template< typename System >
    void operator() ( System &system , state_type &x_tmp , const state_type &x , const state_type *k_vector , double t , const double dt )
    {
        system( x , k_vector[0] , t + m_c * dt );
        algebra<state_type , 1>::foreach( x_tmp , x , m_a_row , k_vector , dt);
    }

private:

    vector< value_type > m_a_row;
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
        m_b = b;
        m_c = c[stage-1];
    }

    template< typename System >
    void operator() ( System &system , state_type &x_tmp , state_type &x , state_type *k_vector , double t , const double dt )
    {
        system( x_tmp , k_vector[stage-1] , t + m_c * dt );
        algebra<state_type , stage>::foreach( x , x , m_b , k_vector , dt);
    }

private:

    vector< value_type > m_b;
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
        m_b = b;
        m_c = c[0];
    }

    template< typename System >
    void do_step( System &system , state_type &x_tmp , state_type &x , state_type k_vector[1] , double t , const double dt )
    {
        system( x , k_vector[0] , t + m_c * dt );
        algebra<state_type , 1>::foreach( x , x , m_b , k_vector , dt);
    }

private:

    vector< value_type > m_b;
    value_type m_c;

};


/* Runge Kutta Stepper - consisting of several stages */

template< typename State , size_t stage_count >
class runge_kutta_stepper
{

	typedef State state_type;
	typedef vector< vector< double > > parameter_array_type2D;
	typedef vector< double > parameter_array_type1D;

    template< class System >
    struct calculate_stage
    {
    	System &system;
    	state_type &x , &x_tmp;
    	state_type *k_vector;
    	const double t;
    	const double dt;
    	const parameter_array_type2D &m_a;
    	const parameter_array_type1D &m_b;
    	const parameter_array_type1D &m_c;

    	calculate_stage( System &_system , state_type &_x , state_type &_x_tmp , state_type *_k_vector ,
    						const double _t , const double _dt ,
    						const parameter_array_type2D &a ,
    						const parameter_array_type1D &b ,
    						const parameter_array_type1D &c )
    	: system( _system ) , x( _x ) , x_tmp( _x_tmp ) , k_vector( _k_vector ) , t( _t ) , dt( _dt ),
    	  m_a( a ) , m_b( b ) , m_c( c )
    	{}

    	template< class Stage >
    	void operator()( Stage &s )
    	{
    		s.init( m_a , m_b , m_c );
    		s( x_tmp , x , k_vector , t , dt );
    	}

    };

public:

    typedef mpl::range_c< size_t , 1 , stage_count > indices;
    typedef typename mpl::copy // create mpl::vector< stepper_stage< 0 >, stepper_stage< 1 > , .. stepper_stage< stage_count-1 > >
    <
      indices ,
      mpl::inserter
      <
        mpl::vector0<> ,
        stepper_stage< State , mpl::_2::value >
      >
    >::type stage_vector;
    typedef stepper_last_stage< State , stage_count > last_stage_type;

    runge_kutta_stepper( const parameter_array_type2D &a ,
							const parameter_array_type1D &b ,
							const parameter_array_type1D &c )
    	: m_a( a ) , m_b( b ) , m_c( c )
    {
    	last_stage.init( b , c );
    }

    template< class System >
    void do_step( System &system , state_type &x , double t , const double dt )
    {
        mpl::for_each< stage_vector >( calculate_stage< System >( system , x , m_x_tmp , m_k_vector , t , dt ,
        													 	   m_a , m_b , m_c ) );
        last_stage.do_step( system , m_x_tmp , x , m_k_vector , t , dt );
    }

private:

	state_type m_k_vector[stage_count];
	state_type m_x_tmp;
	last_stage_type last_stage;
	const parameter_array_type2D &m_a;
	const parameter_array_type1D &m_b;
	const parameter_array_type1D &m_c;
};


#endif /* GENERIC_STEPPER_HPP_ */
