/* Boost odeint/concepts/state_concept.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 This file contains the concepts used in the odeint library

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_CONCEPTS_STATE_CONCEPT_HPP
#define BOOST_NUMERIC_ODEINT_CONCEPTS_STATE_CONCEPT_HPP

#include <boost/concept_check.hpp>

namespace boost {
namespace numeric {
namespace odeint {

    /* Concept StateType */
    template<class X>
    struct Container {

    public:
        typedef typename X::iterator iterator; // requires iterator typedef

        // requires iterator being ForwardIterator
        BOOST_CONCEPT_ASSERT((ForwardIterator<iterator>));

        BOOST_CONCEPT_USAGE(Container)
        {
            same_type(state.begin(), it); //requires begin() method
            same_type(state.end(), it); // requires end() method
        }

    private:

        X state;
        iterator it;

        template<class T>
        void same_type( T const&, T const& );

    };


    /* Concept Resizable */
    template<class X>
    struct Resizer {

    public:

	BOOST_CONCEPT_ASSERT((Container<X>));

	BOOST_CONCEPT_USAGE(Resizer)
	{
	    state.resize(1); // state has resize function
	    size_t n = 2;
	    state.resize(n);
	}

    private:
	X state;

	};

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif
