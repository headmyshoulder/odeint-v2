/* Boost odeint/stepper_rk_generic.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner
 
 This file includes generic Runge Kutta solver
 for ordinary differential equations.

 It solves any ODE dx/dt = f(x,t) using a general Runge Kutta scheme.

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_RK_GENERIC_HPP
#define BOOST_NUMERIC_ODEINT_STEPPER_RK_GENERIC_HPP

#include <vector>
#include <exception>
#include <cmath> // for pow( ) and abs()
#include <limits>

#include <boost/numeric/odeint/container_traits.hpp>


namespace boost {
namespace numeric {
namespace odeint {


    class butcher_tableau_consistency_exception : public std::exception {
        
        virtual const char* what() const throw()
        {
            return "Consistency Check of Butcher Tableau failed!";
        }
    };

    class butcher_tableau_order_condition_exception : public std::exception {

        virtual const char* what() const throw()
        {
            return "Order Condition Check of Butcher Tableau failed!";
        }
    };


    template<
        class Container ,
        class Time = double ,
        class Traits = container_traits< Container >
        >
    class stepper_rk_generic
    {


        // provide basic typedefs
    public:

        typedef unsigned short order_type;
        typedef Container container_type;
        typedef Time time_type;
        typedef Traits traits_type;
        typedef typename traits_type::value_type value_type;
        typedef typename traits_type::iterator iterator;
        typedef typename traits_type::const_iterator const_iterator;





        // private variables
    private:

        typedef std::vector< container_type > container_vector;
        typedef std::vector< iterator > container_iterator_vector;
        typedef std::vector< time_type > parameter_vector;
        typedef std::vector< parameter_vector > parameter_matrix;

        container_vector m_xvec;
        container_iterator_vector m_xiter_vec;
        container_type m_xtmp;
        const parameter_vector m_a;
        const parameter_matrix m_b;
        const parameter_vector m_c;

        order_type m_q;



	// private member functions
    private:

        void reset_iter(typename container_iterator_vector::iterator xiter_iter)
        {
            typename container_vector::iterator x_iter = m_xvec.begin();
            while( x_iter != m_xvec.end() ) {
                (*xiter_iter++) = traits_type::begin(*x_iter++);
            }
        }

        void check_consitency()
        {
            typename parameter_vector::const_iterator a_iter = m_a.begin();
            typename parameter_vector::const_iterator c_iter = m_c.begin();
            typename parameter_matrix::const_iterator b_iter = m_b.begin();
            typename parameter_vector::const_iterator b_iter_iter;

            // check 1: a_i = sum_j b_ij 
            while( a_iter != m_a.end() ) {
                time_type tmp = 0.0;
                b_iter_iter = (*b_iter).begin();
                while( b_iter_iter != (*b_iter).end() ) {
                    tmp += *b_iter_iter++;
                }
                b_iter++;
                if( *a_iter++ != tmp ) 
                    throw butcher_tableau_consistency_exception();
            }

            // check 2: sum_i c_i * (a_i)^(k-1) = 1/k   k = 1..q
            for( unsigned short k = 1; k <= m_q; k++ ) {
                time_type tmp = 0.0;
                a_iter = m_a.begin();
                c_iter = m_c.begin();
                if( k == 1 ) // special treatment for 0^0 = 1
                    tmp += *c_iter++;
                else
                    c_iter++;
                while( a_iter != m_a.end() ) {
                    //std::clog<<pow( *a_iter , k-1 )<< ", ";
                    tmp += (*c_iter++) * pow( *a_iter++ , k-1 );
                }
                //std::clog << tmp << " = " << time_type(1.0)/time_type(k) << "?" << std::endl;
                if( std::abs(time_type(tmp - time_type(1.0)/time_type(k))) > 
                    std::numeric_limits<time_type>::epsilon() ) {
                    //std::clog << std::abs(time_type(tmp - time_type(1.0)/time_type(k))) << std::endl;
                    throw butcher_tableau_order_condition_exception();
                }
            }
        }


    public:



        /* Constructor

           a,b,c are vectors providing the butcher tableau for the Runge Kutta scheme
           
           0     |
           a_1   | b_21
           a_2   | b_31 b_32
           ...   | ...  ...
           a_q-1 | b_q1 b_q2 ... b_qq-1
           -------------------------------
                 | c_1  c_2  ... c_q-1  c_q

          The size of c is determining the order of the scheme q.
          a is of length q-1
          b is of length q-1, b[0] of length 1, b[1] of length 2 and so on
          c is of length q (which defines q).

          There are 2 conditions that these parameters have to fullfill:
          Consitency:

              a_i = sum_j b_ij  for i = 1 ... q-1

          Condition on the order of the method (all error terms dt^k , k<q+1 cancel out):

              sum_i c_i * (a_(i+1))^(k-1) = 1/k   k = 1 ... q

              Note, that a_0 = 1 (implicitely) and 0^0 = 1
              so this means sum_i c_i = 1 at k=1
        */
        stepper_rk_generic( parameter_vector &a,
                            parameter_matrix &b,
                            parameter_vector &c)
            : m_a(a), m_b(b), m_c(c), m_q(c.size())
        {
            m_xvec.resize(m_q);
            m_xiter_vec.resize(m_q);

            check_consitency();
        }

        order_type order() const { return m_q; }


        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                      container_type &x ,
                      container_type &dxdt ,
                      time_type t ,
                      time_type dt )
        {
            using namespace detail::it_algebra;
            typename container_vector::iterator x_iter = m_xvec.begin();
            typename container_iterator_vector::iterator xiter_iter = m_xiter_vec.begin();

            (*x_iter) = dxdt;
            (*xiter_iter++) = traits_type::begin(*x_iter++);

            while( x_iter != m_xvec.end() )
            {
                traits_type::adjust_size(x, (*x_iter));
                (*xiter_iter++) = traits_type::begin(*x_iter++);
            }
            traits_type::adjust_size(x, m_xtmp);
            
            x_iter = m_xvec.begin()+1;
            
            typename parameter_vector::const_iterator a_iter = m_a.begin();
            typename parameter_matrix::const_iterator b_iter = m_b.begin();
            while( x_iter != m_xvec.end() )
            {
                reset_iter(m_xiter_vec.begin());
                scale_sum_generic( traits_type::begin(m_xtmp), traits_type::end(m_xtmp),
                                   (*b_iter).begin(), (*b_iter).end(), dt,
                                   traits_type::begin(x), m_xiter_vec.begin() );
                system( m_xtmp, *x_iter++ , t + dt*(*a_iter++) );
                b_iter++;
            }

            reset_iter(m_xiter_vec.begin());
            typename parameter_vector::const_iterator c_iter = m_c.begin();
            scale_sum_generic( traits_type::begin(x), traits_type::end(x),
                               m_c.begin(), m_c.end(), dt,
                               traits_type::begin(x), m_xiter_vec.begin() );
        }

        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                        container_type &x ,
                        time_type t ,
                        time_type dt )
        {
            traits_type::adjust_size(x, m_xtmp);
            system(x, m_xtmp, t);
            do_step( system, x, m_xtmp, t, dt);
        }

    };



    /* ############################################################
       ############################################################

       C-Array Version of a,b,c handling

       ############################################################
       ############################################################
    */




    template<
        class Container ,
        class Time = double ,
        class Traits = container_traits< Container >
        >
    class stepper_rk_generic_ptr
    {


        // provide basic typedefs
    public:

        typedef const unsigned short order_type;
        typedef Container container_type;
        typedef Time time_type;
        typedef Traits traits_type;
        typedef typename traits_type::value_type value_type;
        typedef typename traits_type::iterator iterator;
        typedef typename traits_type::const_iterator const_iterator;

        // private variables
    private:

        typedef std::vector< container_type > container_vector;
        typedef std::vector< iterator > container_iterator_vector;
        typedef std::vector< time_type > parameter_vector;
        typedef std::vector< parameter_vector > parameter_matrix;

        container_vector m_xvec;
        container_iterator_vector m_xiter_vec;
        container_type m_xtmp;
        const time_type* m_a;
        const time_type* m_b;
        const time_type* m_c;

        order_type m_q;


    private:
        void reset_iter(typename container_iterator_vector::iterator xiter_iter)
        {
            typename container_vector::iterator x_iter = m_xvec.begin();
            while( x_iter != m_xvec.end() ) {
                (*xiter_iter++) = traits_type::begin(*x_iter++);
            }
        }

        void check_consitency()
        {
            const time_type* a_iter = &m_a[0];
            const time_type* b_iter = &m_b[0];
            const time_type* c_iter = &m_c[0];

            unsigned short b_len = 1;
            // check 1: a_i = sum_j b_ij 
            while( a_iter != &m_a[+m_q-1] ) {
                time_type tmp = 0.0;
                const time_type* b_end = b_iter + b_len;
                while( b_iter != b_end ) {
                    tmp += *b_iter++;
                }
                b_len++;
                if( *a_iter++ != tmp ) 
                    throw butcher_tableau_consistency_exception();
            }

            // check 2: sum_i c_i * (a_i)^(k-1) = 1/k   k = 1..q
            for( unsigned short k = 1; k <= m_q; k++ ) {
                time_type tmp = 0.0;
                a_iter = &m_a[0];
                c_iter = &m_c[0];
                if( k == 1 ) // special treatment for 0^0 = 1
                    tmp += *c_iter++;
                else
                    c_iter++;
                while( a_iter != &m_a[m_q-1] ) {
                    tmp += (*c_iter++) * pow( *a_iter++ , k-1 );
                }
                if( std::abs(time_type(tmp - time_type(1.0)/time_type(k))) > 
                    std::numeric_limits<time_type>::epsilon() ) {
                    //std::clog << std::abs(time_type(tmp - time_type(1.0)/time_type(k))) << std::endl;
                    throw butcher_tableau_order_condition_exception();
                }
            }
        }


    public:

        /* Constructor
           
           Same as for the stepper_rk_generic class, but with a, b and c provided 
           as time_type*. The order q of the integration scheme is given separately 
           as parameter. a is a pointer to an array of time_type[] of length q-1,
           b is of length 1 + 2 + ... q-1 = (q-1)(q-2)/2 and ordered as follows:
           b[0] = b_21, b[1] = b_31, b[2] = b_32, b[3] = b_41, b[4] = b42 ...
           c has length q.
           
        */
        stepper_rk_generic_ptr( const time_type* a,
                                const time_type* b,
                                const time_type* c,
                                const unsigned short q)
            : m_a(a), m_b(b), m_c(c), m_q(q)
        {
            m_xvec.resize(m_q);
            m_xiter_vec.resize(m_q);

            check_consitency();
        }

        order_type order() const { return m_q; }

        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                        container_type &x ,
                        container_type &dxdt ,
                        time_type t ,
                        time_type dt )
        {
            using namespace detail::it_algebra;
            typename container_vector::iterator x_iter = m_xvec.begin();
            typename container_iterator_vector::iterator xiter_iter = m_xiter_vec.begin();

            (*x_iter) = dxdt;
            (*xiter_iter++) = traits_type::begin(*x_iter++);

            while( x_iter != m_xvec.end() )
            {
                traits_type::adjust_size(x, (*x_iter));
                (*xiter_iter++) = traits_type::begin(*x_iter++);
            }
            traits_type::adjust_size(x, m_xtmp);
            
            x_iter = m_xvec.begin()+1;
            
            const time_type* a_iter = &m_a[0];
            const time_type* b_iter = &m_b[0];
            unsigned short b_len= 1;
            while( x_iter != m_xvec.end() )
            {
                reset_iter(m_xiter_vec.begin());
                const time_type* b_end = b_iter + b_len;
                scale_sum_generic( traits_type::begin(m_xtmp), traits_type::end(m_xtmp),
                                   b_iter, b_end, dt,
                                   traits_type::begin(x), m_xiter_vec.begin() );
                system( m_xtmp, *x_iter++ , t + dt*(*a_iter++) );
                b_iter = b_end;
                b_len++;
            }

            reset_iter(m_xiter_vec.begin());
            scale_sum_generic( traits_type::begin(x), traits_type::end(x),
                               &m_c[0], &m_c[m_q], dt,
                               traits_type::begin(x), m_xiter_vec.begin() );
        }



        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                        container_type &x ,
                        time_type t ,
                        time_type dt )
        {
            traits_type::adjust_size(x, m_xtmp);
            system(x, m_xtmp, t);
            do_step( system, x, m_xtmp, t, dt);
        }

    };

}
}
}

#endif
