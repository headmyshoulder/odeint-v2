/*
  [auto_generated]
  helpers.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef HELPERS_HPP_INCLUDED
#define HELPERS_HPP_INCLUDED

#include <cstdlib>
#include <algorithm>
#include <utility>

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
std::pair< double , double > tester( Test const& test , size_t n )
{
    timer t;
    double sum = 0.0;
    for( size_t i=0 ; i<n ; ++i )
        sum += test.run();
    return std::make_pair( t.seconds() , sum );
};





#endif // HELPERS_HPP_INCLUDED
