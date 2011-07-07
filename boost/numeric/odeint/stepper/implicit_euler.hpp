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

#include <utility>

#include <boost/ref.hpp>
#include <boost/bind.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

//#include <boost/numeric/odeint/util/size_adjuster.hpp>
//#include <boost/numeric/odeint/util/ublas_resize.hpp>

#include <boost/numeric/odeint/util/ublas_wrapper.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace boost {
namespace numeric {
namespace odeint {








template< class ValueType , class Resizer = initially_resizer >
class implicit_euler
{
	/*void initialize( void )
	{
		m_state_adjuster.register_state( 0 , m_dxdt );
		m_state_adjuster.register_state( 1 , m_x );
		m_state_adjuster.register_state( 2 , m_b );
		m_matrix_adjuster.register_state( 0 , m_jacobi );
		m_pmatrix_adjuster.register_state( 0 , m_pm );

	}

	void copy( const implicit_euler &euler )
	{
		m_dxdt = euler.m_dxdt;
		m_x = euler.m_x;
		m_b = euler.m_b;
		m_jacobi = euler.m_jacobi;
		m_pm = euler.m_pm;
		m_epsilon = euler.m_epsilon;
	}
    */
public:

    typedef ValueType value_type;
    typedef value_type time_type;
    typedef boost::numeric::ublas::vector< value_type > state_type;
    typedef state_wrapper< state_type > wrapped_state_type;
    typedef state_type deriv_type;
    typedef state_wrapper< deriv_type > wrapped_deriv_type;
    typedef boost::numeric::ublas::matrix< value_type > matrix_type;
    typedef state_wrapper< matrix_type > wrapped_matrix_type;
    typedef boost::numeric::ublas::permutation_matrix< size_t > pmatrix_type;
    typedef state_wrapper< pmatrix_type > wrapped_pmatrix_type;
    typedef Resizer resizer_type;
	typedef stepper_tag stepper_category;
	typedef implicit_euler< ValueType , Resizer > stepper_type;

    implicit_euler( const value_type epsilon = 1E-6 )
    : m_epsilon( epsilon ) 
      //m_state_adjuster() , m_matrix_adjuster() , m_pmatrix_adjuster() ,
      //m_dxdt() , m_x() , m_b() ,
      //m_jacobi() , m_pm( 1 )
    {
    	//initialize();
    }
/*
    implicit_euler( const implicit_euler &euler )
    : m_epsilon( 1.0e-6 ) ,
      m_state_adjuster() , m_matrix_adjuster() , m_pmatrix_adjuster() ,
      m_dxdt() , m_x() , m_b() ,
      m_jacobi() , m_pm( 1 )
    {
    	initialize();
    	copy( euler );
    }

    implicit_euler& operator=( const implicit_euler &euler )
    {
    	copy( euler );
    	return *this;
    }
*/
    template< class System >
    void do_step( System system , state_type &x , value_type t , value_type dt )
    {
    	typedef typename boost::unwrap_reference< System >::type system_type;
    	typedef typename boost::unwrap_reference< typename system_type::first_type >::type deriv_func_type;
    	typedef typename boost::unwrap_reference< typename system_type::second_type >::type jacobi_func_type;
    	system_type &sys = system;
    	deriv_func_type &deriv_func = sys.first;
    	jacobi_func_type &jacobi_func = sys.second;

        m_resizer.adjust_size( x , boost::bind( &stepper_type::resize<state_type> , boost::ref( *this ) , _1 ) );

		for( size_t i=0 ; i<x.size() ; ++i )
			m_pm.m_v[i] = i;

        t += dt;

        // apply first Newton step
        deriv_func( x , m_dxdt.m_v , t );

        m_b.m_v = dt * m_dxdt.m_v;

        jacobi_func( x , m_jacobi.m_v  , t );
        m_jacobi.m_v *= dt;
        m_jacobi.m_v -= boost::numeric::ublas::identity_matrix< value_type >( x.size() );

        solve( m_b.m_v , m_jacobi.m_v );

        m_x.m_v = x - m_b.m_v;

        // iterate Newton until some precision is reached
        // ToDo: maybe we should apply only one Newton step -> linear implicit one-step scheme
        while( boost::numeric::ublas::norm_2( m_b.m_v ) > m_epsilon )
        {
            deriv_func( m_x.m_v , m_dxdt.m_v , t );
            m_b.m_v = x - m_x.m_v + dt*m_dxdt.m_v;

            // simplified version, only the first Jacobian is used
//            jacobi( m_x , m_jacobi , t );
//            m_jacobi *= dt;
//            m_jacobi -= boost::numeric::ublas::identity_matrix< value_type >( x.size() );

            solve( m_b.m_v , m_jacobi.m_v );

            m_x.m_v -= m_b.m_v;
        }
        x = m_x.m_v;
    }


    template< class StateIn >
    bool resize( const StateIn &x )
    {
        bool resized = false;
	    resized |= adjust_size_by_resizeability( m_dxdt , x , typename wrapped_deriv_type::is_resizeable() );
        resized |= adjust_size_by_resizeability( m_x , x , typename wrapped_state_type::is_resizeable() );
        resized |= adjust_size_by_resizeability( m_b , x , typename wrapped_deriv_type::is_resizeable() );
        resized |= adjust_size_by_resizeability( m_jacobi , x , typename wrapped_matrix_type::is_resizeable() );
        resized |= adjust_size_by_resizeability( m_pm , x , typename wrapped_pmatrix_type::is_resizeable() );
        return resized;
    }

	template< class StateType >
	void adjust_size( const StateType &x )
	{
		resize( x );
	}


private:

    void solve( state_type &x , matrix_type &m )
    {
        int res = boost::numeric::ublas::lu_factorize( m , m_pm.m_v );
        if( res != 0 ) exit(0);
        boost::numeric::ublas::lu_substitute( m , m_pm.m_v , x );
    }

private:

    value_type m_epsilon;
    resizer_type m_resizer;
    wrapped_deriv_type m_dxdt;
    wrapped_state_type m_x;
    wrapped_deriv_type m_b;
    wrapped_matrix_type m_jacobi;
    wrapped_pmatrix_type m_pm;


};


} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_IMPLICIT_EULER_HPP_INCLUDED
