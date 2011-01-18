/*
 boost header: BOOST_NUMERIC_ODEINT/implicit_euler.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_IMPLICIT_EULER_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_IMPLICIT_EULER_HPP_INCLUDED

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace boost {
namespace numeric {
namespace odeint {


template< class ValueType >
class implicit_euler
{

public:

    typedef ValueType value_type;
    typedef boost::numeric::ublas::vector< value_type > state_type;
    typedef boost::numeric::ublas::matrix< value_type > matrix_type;
    typedef boost::numeric::ublas::permutation_matrix< size_t > pmatrix_type;

    implicit_euler( const value_type epsilon = 1E-6 ) : m_epsilon( epsilon ) , m_pm( 1 )
    { }

    template< class System , class Jacobi >
    void do_step( System system , Jacobi jacobi , state_type &x , value_type t , value_type dt )
    {
		typename boost::unwrap_reference< System >::type &sys = system;
		typename boost::unwrap_reference< Jacobi >::type &jac = jacobi;

        m_dxdt.resize( x.size() );
        m_x.resize( x.size() );
        m_b.resize( x.size() );
        m_jacobi.resize( x.size() , x.size() );
        m_pm = pmatrix_type( x.size() ); // no resize because we also need default filling

        t += dt;

        // apply first Newton step
        sys( x , m_dxdt , t );

        m_b = dt * m_dxdt;

        jac( x , m_jacobi  , t );
        m_jacobi *= dt;
        m_jacobi -= boost::numeric::ublas::identity_matrix< value_type >( x.size() );

        matrix_type jacobi_tmp( m_jacobi );

        solve( m_b , jacobi_tmp );

        m_x = x - m_b;

        // iterate Newton until some precision is reached
        // ToDo: maybe we should apply only one Newton step -> linear implicit one-step scheme
        while( boost::numeric::ublas::norm_2( m_b ) > m_epsilon )
        {
            sys( m_x , m_dxdt , t );
            m_b = x - m_x + dt*m_dxdt;

            /* we use simplified newton where the jacobi matrix is evaluated only once
            jacobi( m_x , m_jacobi , t );
            m_jacobi *= dt;
            m_jacobi -= boost::numeric::ublas::identity_matrix< value_type >( x.size() );
            */
            jacobi_tmp = m_jacobi;

            solve( m_b , jacobi_tmp );

            m_x -= m_b;
        }
        x = m_x;
    }

private:

    void solve( state_type &x , matrix_type &m )
    {
        int res = boost::numeric::ublas::lu_factorize( m , m_pm );
        if( res != 0 ) exit(0);

        boost::numeric::ublas::lu_substitute( m , m_pm , x );
    }

private:
    const value_type m_epsilon;

    state_type m_dxdt;
    state_type m_x;
    state_type m_b;
    matrix_type m_jacobi;
    pmatrix_type m_pm;

};


} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_IMPLICIT_EULER_HPP_INCLUDED
