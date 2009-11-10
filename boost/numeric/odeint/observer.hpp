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

namespace boost {
namespace numeric {
namespace odeint {


    template< class Time , class Container , class System >
    void do_nothing_observer( Time , Container& , System& )
    {
    }
    
   

} // odeint
} // numeric
} // boost


#endif //BOOST_NUMERIC_ODEINT_OBSERVER_HPP_INCLUDED
