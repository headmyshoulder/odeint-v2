/*
 [auto_generated]
 boost/numeric/odeint/integrate/max_step_checker.hpp

 [begin_description]
 Throws exception if too many steps are performed.
 [end_description]

 Copyright 2015 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_MAX_STEP_CHECKER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_INTEGRATE_MAX_STEP_CHECKER_HPP_INCLUDED

#include <stdexcept>

#include <boost/throw_exception.hpp>


namespace boost {
namespace numeric {
namespace odeint {

    /**
     * \brief A class for performing overflow checks on the step count in integrate functions.
     *
     * Provide an instance of this class to integrate functions if you want to throw an exception if
     * too many steps are performed without progress during the integrate routine.
     */
    struct max_step_checker{
        const int m_max_steps;
        int m_steps;

        /**
         * \brief Construct the max_step_checker.
         */
        max_step_checker(const int max_steps = 500)
            : m_max_steps(max_steps)
        {
            reset();
        }

        /**
         * \brief Resets the max_step_checker by setting the internal counter to 0.
         */
        void reset()
        {
            m_steps = 0;
        }

        void operator()(void)
        {
            if( m_steps++ >= m_max_steps )
                BOOST_THROW_EXCEPTION( std::overflow_error( "too many steps!") );
        }
    };

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif