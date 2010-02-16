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


    template< class Time, class Container , class System >
    inline void do_nothing_observer( Time , Container& , System& )
    {
    }
    


    template< class TimeInsertIterator, class StateInsertIterator >
    class state_copy_observer
    {
        
    private:

        TimeInsertIterator m_time_inserter;
        StateInsertIterator m_state_inserter;

    public:

        state_copy_observer( TimeInsertIterator time_inserter ,
			     StateInsertIterator state_inserter ) 
            : m_time_inserter(time_inserter),
	      m_state_inserter(state_inserter)
        {  }

	void reset( void ) { }
        
        template< class Time, class Container, class System >
        void operator () (Time t, Container &state, System &system )
	{
            *m_time_inserter++ = t; // insert time
            *m_state_inserter++ = state; // insert state
        }
    };


} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_OBSERVER_HPP_INCLUDED
