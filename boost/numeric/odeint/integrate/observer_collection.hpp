/*
 * observer_collection.hpp
 *
 *  Created on: Jul 17, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_OBSERVER_COLLECTION_HPP_
#define BOOST_NUMERIC_ODEINT_INTEGRATE_OBSERVER_COLLECTION_HPP_

#include <vector>

#include <boost/function.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template< class State , class Time >
struct observer_collection
{
    typedef boost::function< void( const State& , const Time& ) > observer_type;

    std::vector< observer_type > observers;

    void operator()( const State& x , const Time &t )
    {
        for( size_t i=0 ; i<observers.size() ; ++i )
            observers[i]( x , t );
    }
};

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_INTEGRATE_OBSERVER_COLLECTION_HPP_ */
