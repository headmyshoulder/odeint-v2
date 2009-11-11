/*
 boost header: numeric/odeint/observer.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_OBSERVER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_OBSERVER_HPP_INCLUDED

#include<vector>

namespace boost {
namespace numeric {
namespace odeint {


    template< class Time , class Container , class System >
    void do_nothing_observer( Time , Container& , System& )
    {
    }
    

    template< class InsertIterator, class TimeContainer = std::vector<double> >
    class state_copy_observer {
        
    private:
        TimeContainer *m_times;
        InsertIterator m_state_inserter;
        typename TimeContainer::iterator m_time_iter;

        typedef typename TimeContainer::value_type time_type;

    public:
        state_copy_observer( TimeContainer &times, InsertIterator state_inserter ) 
            : m_times(&times), m_state_inserter(state_inserter), m_time_iter(m_times->begin()) 
        {  }
        
        template< class Container, class System >
        void operator () (time_type t, Container &state, System &system ) {
            if( t >= *m_time_iter ) { // we've reached the next time point
                *m_state_inserter++ = state; // insert the state
                m_time_iter++; // next time point
            }
        }
            
    };
   

} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_OBSERVER_HPP_INCLUDED
