

#ifndef BOOST_NUMERIC_ODEINT_ERROR_CHECKER_STANDARD_HPP
#define BOOST_NUMERIC_ODEINT_ERROR_CHECKER_STANDARD_HPP

#include <cmath>

namespace boost {
namespace numeric {
namespace odeint {

    template< class Container, class Time >
    class error_checker_standard {


    public:
        typedef Container container_type;
        typedef Time time_type;
        typedef typename container_type::value_type value_type;
        typedef typename container_type::iterator iterator;

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
            iterator x_iter = x.begin();
            const iterator x_end = x.end();
            iterator dxdt_iter = dxdt.begin();
            iterator scale_iter = scale.begin();
            while( x_iter != x_end ) 
            {
                *scale_iter++ = m_eps_abs + 
                    m_eps_rel * (m_a_x * std::abs(*x_iter++) + 
                                 m_a_dxdt * dt * std::abs(*dxdt_iter++));
            }
        }

        time_type get_max_error_ratio( container_type &x_err, container_type &scale)
        {
            iterator x_err_iter = x_err.begin();
            const iterator x_err_end = x_err.end();
            iterator scale_iter = scale.begin();
            time_type max_rel_error = static_cast<time_type>(0.0);
            while( x_err_iter != x_err_end ) 
            {
                max_rel_error = std::max(
                        static_cast<time_type>(std::abs(*x_err_iter++)/(*scale_iter++)) , 
                        max_rel_error);
            }            
            return max_rel_error;
        }

        const time_type get_epsilon() { return std::max(m_eps_rel, m_eps_abs); }


    };

}
}
}

#endif
