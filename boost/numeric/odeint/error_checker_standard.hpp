/* Boost odeint/error_checker_standard.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file includes the standard error checker to be used with
 controlled steppers. It's purpose is to provide a method
 that calculates the error ration of a given error-estimation
 with respect to some user defined tolerance.

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_ERROR_CHECKER_STANDARD_HPP
#define BOOST_NUMERIC_ODEINT_ERROR_CHECKER_STANDARD_HPP

#include <cmath>

namespace boost {
namespace numeric {
namespace odeint {

    template< class Container, 
              class Time, 
              class Traits = container_traits< Container > >
    class error_checker_standard {


    public:
        typedef Time time_type;
        typedef Traits traits_type;
        typedef typename traits_type::container_type container_type;
        typedef typename traits_type::value_type value_type;
        typedef typename traits_type::iterator iterator;
        typedef typename traits_type::const_iterator const_iterator;


    private:
       	time_type m_eps_abs;
        time_type m_eps_rel;
        time_type m_a_x;
        time_type m_a_dxdt;

    public:
        // constructor
	error_checker_standard( 
                time_type abs_err, time_type rel_err, 
                time_type factor_x, time_type factor_dxdt )
            : m_eps_abs( abs_err ), m_eps_rel( rel_err ),
              m_a_x( factor_x ), m_a_dxdt( factor_dxdt )
        { }

        void fill_scale( 
                container_type &x, 
                container_type &dxdt, 
                time_type dt, 
                container_type &scale )
        {
            detail::it_algebra::weighted_error( traits_type::begin(scale),
                                                traits_type::end(scale),
                                                traits_type::begin(x),
                                                traits_type::end(dxdt),
                                                m_eps_abs,
                                                m_eps_rel,
                                                m_a_x,
                                                m_a_x*dt );
        }

        time_type get_max_error_ratio( container_type &x_err, container_type &scale)
        {
            return detail::it_algebra::max_ratio( traits_type::begin(x_err),
                                                  traits_type::end(x_err),
                                                  traits_type::begin(scale),
                                                  static_cast<time_type>(0.0) );
        }

        const time_type get_epsilon() { return std::max(m_eps_rel, m_eps_abs); }


    };

}
}
}

#endif
