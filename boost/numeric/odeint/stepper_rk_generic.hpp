/* Boost odeint/stepper_rk_generic.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner
 
 This file includes generic Runge Kutta solver
 for ordinary differential equations.

 It solves any ODE dx/dt = f(x,t).

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_RK_GENERIC_HPP
#define BOOST_NUMERIC_ODEINT_STEPPER_RK_GENERIC_HPP

#include <vector>
#include <boost/concept_check.hpp>

#include <boost/numeric/odeint/concepts/state_concept.hpp>
#include <boost/numeric/odeint/resizer.hpp>

#include <iostream>

namespace boost {
namespace numeric {
namespace odeint {

    template<
        class Container ,
        class Time = double ,
        class Resizer = resizer< Container >
        >
    class stepper_rk_generic {

        // provide basic typedefs
    public:

        typedef Container container_type;
        typedef Resizer resizer_type;
        typedef Time time_type;
        typedef const unsigned short order_type;
        typedef typename container_type::value_type value_type;
        typedef typename container_type::iterator iterator;


        // check the concept of the ContainerType
    private:

        BOOST_CLASS_REQUIRE( container_type ,
                             boost::numeric::odeint, Container );

        // private variables
    private:
        typedef std::vector< container_type > container_vector;
        typedef std::vector< iterator > container_iterator_vector;

        container_vector m_xvec;
        container_iterator_vector m_xiter_vec;
        container_type m_xtmp;
        const std::vector< time_type > m_a;
        const std::vector< std::vector<time_type> > m_b;
        const std::vector< time_type > m_c;

        order_type m_q;

        resizer_type m_resizer;

    public:

        stepper_rk_generic( std::vector<time_type> &a,
                            std::vector< std::vector<time_type> > &b,
                            std::vector< time_type > &c)
            : m_a(a), m_b(b), m_c(c), m_q(c.size())
        {
            m_xvec.resize(m_q);
            m_xiter_vec.resize(m_q);
        }

        order_type order() const { return m_q; }

        template< class DynamicalSystem >
        void next_step( DynamicalSystem &system ,
                        container_type &x ,
                        container_type &dxdt ,
                        time_type t ,
                        time_type dt )
        {
            using namespace detail::it_algebra;
            typename container_vector::iterator x_iter = m_xvec.begin();
            typename container_iterator_vector::iterator xiter_iter = m_xiter_vec.begin();
            while( x_iter != m_xvec.end() ) {
                m_resizer.adjust_size(x, (*x_iter));
                (*xiter_iter++) = (*x_iter).begin();
                x_iter++;
            }
            m_resizer.adjust_size(x, m_xtmp);
            
            x_iter = m_xvec.begin(); 
            (*x_iter++) = dxdt;
            
            typename std::vector< time_type >::const_iterator a_iter = m_a.begin();
            typename std::vector< std::vector<time_type> >::const_iterator b_iter = m_b.begin();
            while( x_iter != m_xvec.end() ) {
                reset_iter(m_xiter_vec.begin());
                scale_sum_generic( m_xtmp.begin(), m_xtmp.end(),
                                   (*b_iter).begin(), (*b_iter).end(), dt,
                                   x.begin(), m_xiter_vec.begin() );
                system( m_xtmp, *x_iter , t + dt*(*a_iter) );
                x_iter++;
                a_iter++;
                b_iter++;
            }

            reset_iter(m_xiter_vec.begin());
            typename std::vector< time_type >::const_iterator c_iter = m_c.begin();
            scale_sum_generic( x.begin(), x.end(),
                               m_c.begin(), m_c.end(), dt,
                               x.begin(), m_xiter_vec.begin() );
        }

        void reset_iter(typename container_iterator_vector::iterator xiter_iter)
        {
            typename container_vector::iterator x_iter = m_xvec.begin();
            while( x_iter != m_xvec.end() ) {
                (*xiter_iter++) = (*x_iter++).begin();
            }
        }

        template< class DynamicalSystem >
        void next_step( DynamicalSystem &system ,
                        container_type &x ,
                        time_type t ,
                        time_type dt )
        {
            m_resizer.adjust_size(x, m_xtmp);
            system(x, m_xtmp, t);
            next_step( system, x, m_xtmp, t, dt);
        }

    };

}
}
}

#endif
