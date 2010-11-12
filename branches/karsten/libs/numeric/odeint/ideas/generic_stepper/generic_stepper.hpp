/*
 * generic_stepper.hpp
 *
 *  Created on: Nov 12, 2010
 *      Author: mario
 */

#ifndef GENERIC_STEPPER_HPP_
#define GENERIC_STEPPER_HPP_


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
        m_a_row( a[stage-1] );
        m_c( c[stage] );
    }

    operator() ( state_type &x_tmp , const state_type &x , const state_type *k_vector , double t , const double dt )
    {
        system( x_tmp , k_vector[stage] , t + m_c * dt );
        algebra::foreach< stage >( x_tmp , x , m_a_row , k_vector , dt);
    }

private:

    vector< value_type > m_a_row;
    value_type m_c;

};


template<>
class stepper_stage< 0 >
{

public:

    typedef double value_type;
    typedef State state_type;
    typedef vector< vector< double > > parameter_array_type2D;
    typedef vector< double > parameter_array_type1D;


    void init( parameter_array_type2D &a , parameter_array_type1D &c )
    {
        m_c( c[0] );
    }

    operator() ( state_type &x_tmp , const state_type &x , const state_type *k_vector , double t , const double dt )
    {
        system( x , k_vector[0] , t + c * dt );
        algebra::foreach< stage >( x_tmp , x , m_a_row , k_vector , dt);
    }

private:

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


    void init( parameter_array_type2D &a , parameter_array_type1D &c )
    {
        m_a_row( a[stage-1] );
        m_c( c[stage] );
    }

    operator() ( state_type &x_tmp , const state_type &x , const state_type *k_vector , double t , const double dt )
    {
        system( x_tmp , k_vector[stage] , t + m_c * dt );
        algebra::foreach< stage >( x , x , m_a_row , k_vector , dt);
    }

private:

    vector< value_type > m_a_row;
    value_type m_c;

};


/* Runge Kutta Stepper - consisting of several stages */

template< size_t stage_count >
class runge_kutta_stepper
{
    struct initialize
    {
        typedef double value_type;
        typedef vector< vector< double > > parameter_array_type2D;
        typedef vector< double > parameter_array_type1D;

        parameter_array_type2D &a;
        parameter_array_type1D &c;

        initialize( parameter_array_type2D &_a , parameter_array_type1D &_c ) : a( _a ) , c( _c ) { }


        template< typename Stage >
        operator() ( Stage &s ) { s.init( a , c ); }
    };

    template< typename System >
    struct calculate_stage
    {
        calculate_stage() {}
    };

public:

    typedef mpl::range_c< size_t , 0 , stage_count-1 > indices;
    typedef mpl::copy // create mpl::vector< stepper_stage< 0 >, stepper_stage< 1 > , .. stepper_stage< stage_count-1 > >
    <
      indices ,
      mpl::inserter
      <
        mpl::vector0<> ,
        mpl::insert_range
        <
          mpl::_1 ,
          mpl::end< mpl::_1 > ,
          stage< mpl::_2 >
        >
      >
    >::type stages;
    typedef stepper_last_stage< stage_count > last_stage;

    runge_kutta_stepper( parameter_array_type2D &a , parameter_array_type1D &c ) : last_stage( a , c )
    {
        mpl::for_each< stages >( initialize( a , c ) );
    }

    void do_step()
    {
        mpl::for_each< stages >( calculate_stage< System >( system , x , m_x_tmp , m_k_vector , t , dt ) );
    }

};


#endif /* GENERIC_STEPPER_HPP_ */
