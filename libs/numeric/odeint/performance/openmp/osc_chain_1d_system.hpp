/* Boost libs/numeric/odeint/performance/openmp/osc_chain_1d_system.hpp

 Copyright 2009-2012 Karsten Ahnert
 Copyright 2009-2012 Mario Mulansky

 stronlgy nonlinear hamiltonian lattice

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <vector>
#include <cmath>
#include <iostream>

#include <omp.h>

#include <boost/math/special_functions/sign.hpp>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>

typedef std::vector< double > dvec;

namespace checked_math {
    inline double pow( double x , double y )
    {
        if( x==0.0 )
            // 0**y = 0, don't care for y = 0 or NaN
            return 0.0;
        using std::pow;
        using std::abs;
        return pow( abs(x) , y );
    }
}

double signed_pow( double x , double k )
{
    using boost::math::sign;
    return checked_math::pow( x , k ) * sign(x);
}

struct osc_chain {

    const double m_kap, m_lam, m_beta;

    osc_chain( const double kap , const double lam , const double beta )
        : m_kap( kap ) , m_lam( lam ) , m_beta( beta )
    { }

    // Simple case with openmp_range_algebra
    void operator()( const std::vector<double> &q ,
                           std::vector<double> &dpdt ) const
    {
        const size_t N = q.size();
        double coupling_lr = 0;
        size_t last_i = N;
#       pragma omp parallel for firstprivate(coupling_lr, last_i) schedule(runtime)
        for(size_t i = 0 ; i < N ; ++i)
        {
            if(i > 0 && i != last_i + 1)
                coupling_lr = signed_pow( q[i-1]-q[i] , m_lam-1 );
            dpdt[i] = -signed_pow( q[i] , m_kap-1 ) + coupling_lr;
            coupling_lr = signed_pow( q[i] - q[i+1] , m_lam-1 );
            dpdt[i] -= coupling_lr;
            last_i = i;
        }
    }

    // Split case with openmp_algebra
    void operator()( const boost::numeric::odeint::openmp_state<double> &q ,
                           boost::numeric::odeint::openmp_state<double> &dpdt ) const
    {
        const size_t N = q.size();
#       pragma omp parallel for schedule(runtime)
        for(size_t i = 0 ; i < N ; ++i)
        {
            const double q_left =      i == 0 ? 0 : q[i-1].back();
            const double q_right = i + 1 == N ? 0 : q[i+1].front();

            const std::vector<double> &_q = q[i];
            std::vector<double> &_dpdt = dpdt[i];

            const size_t M = q[i].size();
            double coupling_lr = signed_pow( q_left - _q[0] , m_lam-1 );
            for(size_t i = 0 ; i < M-1 ; ++i)
            {
                _dpdt[i] = -signed_pow( _q[i] , m_kap-1 ) + coupling_lr;
                coupling_lr = signed_pow( _q[i] - _q[i+1] , m_lam-1 );
                _dpdt[i] -= coupling_lr;
            }
            _dpdt[N-1] = -signed_pow( _q[N-1] , m_kap-1 ) + coupling_lr;
            _dpdt[N-1] -= signed_pow( _q[N-1] - q_right , m_lam-1 );
        }
    }

};

#endif
