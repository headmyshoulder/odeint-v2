/*
 * const_step_iterator.cpp
 *
 *  Created on: Jun 26, 2012
 *      Author: karsten
 */

#include <iostream>
#include <iterator>
#include <utility>
#include <algorithm>
#include <array>

#include <boost/range/algorithm.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/numeric.hpp>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/iterator/const_step_iterator.hpp>

#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

struct lorenz
{
    template< class State , class Deriv >
    void operator()( const State &x , Deriv &dxdt , double t ) const
    {
        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = -b * x[2] + x[0] * x[1];
    }
};



int main( int argc , char **argv )
{
    typedef std::array< double , 3 > state_type;

    // std::for_each
    {
        runge_kutta4< state_type > stepper;
        state_type x = {{ 10.0 , 10.0 , 10.0 }};
        std::for_each( make_const_step_iterator( stepper , lorenz() , x , 0.0 , 0.01 ) ,
                       make_const_step_iterator( stepper , lorenz() , x , 1.0 , 0.01 ) ,
                       []( const state_type &x ) {
                           std::cout << x[0] << tab << x[1] << tab << x[2] << "\n"; } );
    }

    // std::copy_if
    {
        std::vector< state_type > res;
        runge_kutta4< state_type > stepper;
        state_type x = {{ 10.0 , 10.0 , 10.0 }};
        std::copy_if( make_const_step_iterator( stepper , lorenz() , x , 0.0 , 0.01 ) ,
                      make_const_step_iterator( stepper , lorenz() , x , 1.0 , 0.01 ) ,
                      std::back_inserter( res ) ,
                      []( const state_type &x ) {
                          return ( x[0] > 0.0 ) ? true : false; } );
        for( size_t i=0 ; i<res.size() ; ++i )
            cout << res[i][0] << tab << res[i][1] << tab << res[i][2] << "\n";
    }

    // std::accumulate
    {
        runge_kutta4< state_type > stepper;
        state_type x = {{ 10.0 , 10.0 , 10.0 }};
        double res = std::accumulate( make_const_step_iterator( stepper , lorenz() , x , 0.0 , 0.01 ) ,
                                      make_const_step_iterator( stepper , lorenz() , x , 1.0 , 0.01 ) ,
                                      0.0 ,
                                      []( double sum , const state_type &x ) {
                                          return sum + x[0]; } );
        cout << res << endl;
    }


    // std::transform
    {
        runge_kutta4< state_type > stepper;
        state_type x = {{ 10.0 , 10.0 , 10.0 }};
        vector< double > weights;
        std::transform( make_const_step_iterator( stepper , lorenz() , x , 0.0 , 0.01 ) ,
                        make_const_step_iterator( stepper , lorenz() , x , 1.0 , 0.01 ) ,
                        back_inserter( weights ) ,
                        []( const state_type &x ) {
                            return sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] ); } );
        for( size_t i=0 ; i<weights.size() ; ++i )
            cout << weights[i] << "\n";
    }





    /*
     * Boost.Range versions
     */


    // boost::range::for_each
    {
        runge_kutta4< state_type > stepper;
        state_type x = {{ 10.0 , 10.0 , 10.0 }};
        boost::range::for_each( make_const_step_range( stepper , lorenz() , x , 0.0 , 1.0 , 0.01 ) ,
                                []( const state_type &x ) {
                                    std::cout << x[0] << tab << x[1] << tab << x[2] << "\n"; } );
    }

    // boost::range::copy with filtered adaptor (simulating std::copy_if)
    {
        std::vector< state_type > res;
        runge_kutta4< state_type > stepper;
        state_type x = {{ 10.0 , 10.0 , 10.0 }};
        boost::range::copy( make_const_step_range( stepper , lorenz() , x , 0.0 , 1.0 , 0.01 ) |
                            boost::adaptors::filtered( [] ( const state_type &x ) { return ( x[0] > 0.0 ); } ) ,
                            std::back_inserter( res ) );
        for( size_t i=0 ; i<res.size() ; ++i )
            cout << res[i][0] << tab << res[i][1] << tab << res[i][2] << "\n";
    }

    // boost::range::accumulate
    {
        runge_kutta4< state_type > stepper;
        state_type x = {{ 10.0 , 10.0 , 10.0 }};
        double res = boost::accumulate( make_const_step_range( stepper , lorenz() , x , 0.0 , 1.0 , 0.01 ) , 0.0 ,
                                        []( double sum , const state_type &x ) {
                                            return sum + x[0]; } );
        cout << res << endl;
    }

    //  boost::range::transform
    {
        runge_kutta4< state_type > stepper;
        state_type x = {{ 10.0 , 10.0 , 10.0 }};
        vector< double > weights;
        boost::transform( make_const_step_range( stepper , lorenz() , x , 0.0 , 1.0 , 0.01 ) , back_inserter( weights ) ,
                          []( const state_type &x ) {
                              return sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] ); } );
        for( size_t i=0 ; i<weights.size() ; ++i )
            cout << weights[i] << "\n";
    }



    return 0;
}
