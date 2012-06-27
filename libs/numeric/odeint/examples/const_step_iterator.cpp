/*
 * const_step_ode_iterator.cpp
 *
 *  Created on: Jun 26, 2012
 *      Author: karsten
 */

#include <iostream>
#include <utility>
#include <algorithm>

#include <boost/array.hpp>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/iterator/const_step_iterator.hpp>


using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

struct lorenz
{
    template< class State , class Deriv >
    void operator()( const State &x_ , Deriv &dxdt_ , double t ) const
    {
        typename boost::range_iterator< const State >::type x = boost::begin( x_ );
        typename boost::range_iterator< Deriv >::type dxdt = boost::begin( dxdt_ );

        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = -b * x[2] + x[0] * x[1];
    }
};


struct printer
{
    std::ostream &m_out ;
    printer( std::ostream &out = std::cout ) : m_out( out ) { }

    template< class State >
    void operator()( const State &s )
    {
        m_out << s[0] << "\t" << s[1] << "\t" << s[2] << "\n";
    }
};

int main( int argc , char **argv )
{
    typedef boost::array< double , 3 > state_type;
    runge_kutta4< state_type > stepper;
    state_type x = {{ 10.0 , 10.0 , 10.0 }};
    std::for_each( make_const_step_iterator( stepper , lorenz() , x , 0.0 , 0.01 ) ,
                   make_const_step_iterator( stepper , lorenz() , x , 10.0 , 0.01 ) , printer() );

    return 0;
}
