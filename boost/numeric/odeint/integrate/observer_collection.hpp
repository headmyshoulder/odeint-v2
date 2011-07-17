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
class observer_collection
{
public:

    typedef boost::function< void( const State& , const Time& ) > observer_type;
    typedef std::vector< observer_type > observer_collection_type;

    void operator()( const State& x , const Time &t )
    {
        for( size_t i=0 ; i<m_observers.size() ; ++i )
            m_observers[i]( x , t );
    }

    observer_collection_type& observers( void ) { return m_observers; }
    const observer_collection_type& observers( void ) const { return m_observers; }

private:

    observer_collection_type m_observers;
};

} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_INTEGRATE_OBSERVER_COLLECTION_HPP_ */
