/*
 [auto_generated]
 const_step_range.cpp

 [begin_description]
 tba.
 [end_description]

 Copyright 2009-2012 Karsten Ahnert
 Copyright 2009-2012 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <phase_lattice.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/range/const_step_range.hpp>
#include <boost/numeric/odeint/iterator/const_step_iterator.hpp>
#include <boost/numeric/odeint/iterator/n_step_iterator.hpp>

#include <boost/range/algorithm/for_each.hpp>
#include <boost/range.hpp>

#include <iostream>
#include <algorithm>
#include <fstream>

namespace odeint = boost::numeric::odeint;

const double t_start = 0.0;
const double t_end = 1.0;
const double dt = 0.01;

template< typename State >
void ic( State &x )
{
    srand48( 123 );
    std::generate( x.begin() , x.end() , []() -> double { return 2.0 * 3.1415927 * drand48(); } );
}

struct observer
{
    std::ostream &m_out;
    observer( std::ostream &out ) : m_out( out ) {}
    
    template< typename State >
    void operator()( State const& x ) const
    {
        for( size_t i=0 ; i<x.size() ; ++i )
            m_out << x[i] << " ";
        m_out << "\n";
    }
    
    template< typename State , typename Time >
    void operator()( State const& x , Time const& t ) const
    {
        (*this)( x );
    }
};



template< size_t N >
struct test_range
{
    static void run( void )
    {
        typedef phase_lattice< N > system_type;
        typedef typename system_type::state_type state_type;

        state_type x;
        ic( x );
        std::ofstream fout( "dat/range.dat" );
        boost::for_each( odeint::make_const_step_range(
            odeint::runge_kutta4< state_type > {} ,
            x ,
            phase_lattice< N > {} , 
            t_start , t_end , dt ) , observer { fout } );
    }
};

template< size_t N >
struct test_iterator
{
    static void run( void )
    {
        typedef phase_lattice< N > system_type;
        typedef typename system_type::state_type state_type;

        state_type x;
        ic( x );
        std::ofstream fout( "dat/iterator.dat" );
        std::for_each(
            odeint::make_const_step_iterator_begin(
                odeint::runge_kutta4< state_type > {} ,
                phase_lattice< N > {} , 
                x , t_start , t_end , dt ) ,
            odeint::make_const_step_iterator_end(
                odeint::runge_kutta4< state_type > {} ,
                phase_lattice< N > {} ,
                x ) ,                
            observer { fout } );
    }

};

template< size_t N >
struct test_iterator_range
{
    static void run( void )
    {
        typedef phase_lattice< N > system_type;
        typedef typename system_type::state_type state_type;

        state_type x;
        ic( x );
        std::ofstream fout( "dat/iterator_range.dat" );
        boost::for_each(
            odeint::make_const_step_range(
                odeint::runge_kutta4< state_type > {} ,
                phase_lattice< N > {} , 
                x , t_start , t_end , dt ) ,
            observer { fout } );
    }
};

template< size_t N >
struct test_n_step_iterator_range
{
    static void run( void )
    {
        typedef phase_lattice< N > system_type;
        typedef typename system_type::state_type state_type;

        state_type x;
        ic( x );
        std::ofstream fout( "dat/n_step_iterator_range.dat" );
        boost::for_each(
            odeint::make_n_step_range(
                odeint::runge_kutta4< state_type > {} ,
                phase_lattice< N > {} , 
                x , t_start , dt , 100 ) ,
            observer { fout } );
    }
};


template< size_t N >
struct test_raw
{
    static void run( void )
    {
        typedef phase_lattice< N > system_type;
        typedef typename system_type::state_type state_type;

        state_type x;
        ic( x );
        std::ofstream fout( "dat/raw.dat" );
        observer obs( fout );
        double t = t_start;
        odeint::runge_kutta4< state_type > stepper;
        obs( x );
        while( t < t_end )
        {
            stepper.do_step( phase_lattice< N > {} , x , t , dt );
            t += dt;
            obs( x );
        }
    }
};

template< size_t N >
struct test_integrate
{
    static void run( void )
    {
        typedef phase_lattice< N > system_type;
        typedef typename system_type::state_type state_type;

        state_type x;
        ic( x );
        std::ofstream fout( "dat/integrate.dat" );
        odeint::integrate_const( odeint::runge_kutta4< state_type > {} , phase_lattice< N > {} ,
                                 x , t_start , t_end , dt , observer { fout } );
    }

};

template< size_t N >
struct test
{
    static void run( void )
    {
        test_range< N >::run();
        test_iterator< N >::run();
        test_iterator_range< N >::run();
        test_raw< N >::run();
        test_integrate< N >::run();
        test_n_step_iterator_range< N >::run();
    }
};


int main( int argc , char *argv[] )
{
    test< 8 >::run();
   
    return 0;
}
