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
#include <chrono>

namespace odeint = boost::numeric::odeint;

// Only for arrays

#include <chrono>


class timer
{
    typedef std::chrono::high_resolution_clock clock_type;
    clock_type::time_point m_start_time;

    template< class T >
    static inline double get_seconds( T t )
    {
        return double( std::chrono::duration_cast< std::chrono::nanoseconds >( t ).count() ) * 1.0e-6;
    }

public:

    timer( void ) : m_start_time( clock_type::now() ) { }

    inline double seconds( void ) const
    {
        return get_seconds( clock_type::now() - m_start_time );
    }

    inline void restart( void )
    {
        m_start_time = clock_type::now();
    }
};



template< typename State >
void ic( State &x )
{
    srand48( 123 );
    std::generate( x.begin() , x.end() , []() -> double { return 2.0 * 3.1415927 * drand48(); } );
}

struct empty_observer
{
    template< typename State >
    void operator()( State const& x ) const
    {
    }
    
    template< typename State , typename Time >
    void operator()( State const& x , Time const& t ) const
    {
    }
};

template< typename Test >
std::pair< double , double > tester( size_t n )
{
    timer t;
    double sum = 0.0;
    for( size_t i=0 ; i<n ; ++i )
        sum += Test::run();
    return std::make_pair( t.seconds() , sum );
};


template< size_t N >
struct test_range
{
    static double run( void )
    {
        typedef phase_lattice< N > system_type;
        typedef typename system_type::state_type state_type;

        state_type x;
        ic( x );
        boost::for_each( odeint::make_const_step_range(
            odeint::runge_kutta4< state_type > {} ,
            x ,
            phase_lattice< N > {} , 
            0.0 , 100.0 , 0.01 ) , empty_observer {} );
        return x[0];
    }
};

template< size_t N >
struct test_iterator
{
    static double run( void )
    {
        typedef phase_lattice< N > system_type;
        typedef typename system_type::state_type state_type;

        state_type x;
        ic( x );
        std::for_each(
            odeint::make_const_step_iterator_begin(
                odeint::runge_kutta4< state_type > {} ,
                phase_lattice< N > {} , 
                x , 0.0 , 100.0 , 0.01 ) ,
            odeint::make_const_step_iterator_end(
                odeint::runge_kutta4< state_type > {} ,
                phase_lattice< N > {} ,
                x ) ,                
            empty_observer {} );
        return x[0];
    }

};

template< size_t N >
struct test_iterator_range
{
    static double run( void )
    {
        typedef phase_lattice< N > system_type;
        typedef typename system_type::state_type state_type;

        state_type x;
        ic( x );
        boost::for_each(
            odeint::make_const_step_range(
                odeint::runge_kutta4< state_type > {} ,
                phase_lattice< N > {} , 
                x , 0.0 , 100.0 , 0.01 ) ,
            empty_observer {} );
        return x[0];
    }
};

template< size_t N >
struct test_n_step_iterator_range
{
    static double run( void )
    {
        typedef phase_lattice< N > system_type;
        typedef typename system_type::state_type state_type;

        state_type x;
        ic( x );
        boost::for_each(
            odeint::make_n_step_range(
                odeint::runge_kutta4< state_type > {} ,
                phase_lattice< N > {} , 
                x , 0.0 , 0.01 , 10000 ) ,
            empty_observer {} );
        return x[0];
    }
};


template< size_t N >
struct test_raw
{
    static double run( void )
    {
        typedef phase_lattice< N > system_type;
        typedef typename system_type::state_type state_type;

        state_type x;
        ic( x );
        double t = 0.0;
        odeint::runge_kutta4< state_type > stepper;
        while( t < 100.0 )
        {
            stepper.do_step( phase_lattice< N > {} , x , t , 0.01 );
            t += 0.01;
        }
        return x[0];
    }
};

template< size_t N >
struct test_integrate
{
    static double run( void )
    {
        typedef phase_lattice< N > system_type;
        typedef typename system_type::state_type state_type;

        state_type x;
        ic( x );
        odeint::integrate_const( odeint::runge_kutta4< state_type > {} , phase_lattice< N > {} ,
                                 x , 0.0 , 100.0 , 0.01 , empty_observer {} );
        return x[0];
    }

};

template< size_t N >
struct test
{
    static void run( void )
    {
        auto res0 = tester< test_range< N > >( 10 );
        auto res1 = tester< test_iterator< N > >( 10 );
        auto res2 = tester< test_iterator_range< N > >( 10 );
        auto res3 = tester< test_raw< N > >( 10 );
        auto res4 = tester< test_integrate< N > >( 10 );
        auto res5 = tester< test_n_step_iterator_range< N > >( 10 );
        std::cout << N << " "
                  << res0.first << " " << res1.first << " " << res2.first << " " << res3.first << " " << res4.first << " " << res5.first << " "
                  << res0.second << " " << res1.second << " " << res2.second << " " << res3.second << " " << res4.second << " " << res5.second
                  << std::endl;
    }
};


int main( int argc , char *argv[] )
{
    std::cout.precision( 14 );
    test< 8 >::run();
    test< 16 >::run();
    test< 32 >::run();
    test< 64 >::run();
    test< 128 >::run();
    test< 256 >::run();
    test< 512 >::run();
    test< 1024 >::run();
    test< 2048 >::run();
    test< 4096 >::run();
    
    return 0;
}
