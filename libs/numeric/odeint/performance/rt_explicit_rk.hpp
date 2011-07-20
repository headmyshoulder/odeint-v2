/*
 * rt_explicit_rk.hpp
 *
 *  Created on: May 11, 2011
 *      Author: mario
 */
#ifndef RT_EXPLICIT_RK_HPP_
#define RT_EXPLICIT_RK_HPP_

#include <vector>

#include "rt_algebra.hpp"

using namespace std;

template< class StateType >
class rt_explicit_rk
{
public:
    typedef StateType state_type;
    typedef vector< vector< double > > coeff_a_type;
    typedef vector< double > coeff_b_type;
    typedef vector< double > coeff_c_type;

    rt_explicit_rk( size_t stage_count ) : m_s( stage_count )
    {
        m_F = new state_type[ m_s ];
    }

    rt_explicit_rk( size_t stage_count , coeff_a_type &a , coeff_b_type &b , coeff_c_type &c )
        : m_s( stage_count ) , m_a( a ) , m_b( b ) , m_c( c )
    {
        m_F = new state_type[ m_s ];
    }

    ~rt_explicit_rk()
    {
        delete[] m_F;
    }

    void set_params( coeff_a_type &a , coeff_b_type &b , coeff_c_type &c )
    {
        m_a = a;
        m_b = b;
        m_c = c;
    }

    template< class System >
    void do_step( System sys , state_type &x , const double t , const double dt )
    {
        // first stage separately
        sys( x , m_F[0] , t + m_c[0]*t );
        if( m_s == 1 )
            rt_algebra::foreach( x , x , m_b , m_F , dt );
        else
            rt_algebra::foreach( m_x_tmp , x , m_a[0] , m_F , dt );

        for( size_t stage = 2 ; stage <= m_s ; ++stage )
        {
            sys( m_x_tmp , m_F[stage-1] , t + m_c[stage-1]*dt );
            if( stage == m_s )
                rt_algebra::foreach( x , x , m_b , m_F , dt );
            else
                rt_algebra::foreach( m_x_tmp , x , m_a[stage-1] , m_F , dt );
        }
    }


private:
    size_t m_s;
    vector< vector< double > > m_a;
    vector< double > m_b;
    vector< double > m_c;

    state_type m_x_tmp;
    state_type *m_F;
};

#endif /* RT_EXPLICIT_RK_HPP_ */
