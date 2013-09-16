/* Boost libs/numeric/odeint/performance/openmp/osc_chain_1d.cpp

 Copyright 2009-2012 Karsten Ahnert
 Copyright 2009-2012 Mario Mulansky

 stronlgy nonlinear hamiltonian lattice in 2d

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifdef _OPENMP

#include <iostream>
#include <vector>
#include <random>

#include <omp.h>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>

#include <boost/timer/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/min.hpp>

#include "osc_chain_1d_system.hpp"

using namespace std;
using namespace boost::numeric::odeint;
using namespace boost::accumulators;

using boost::timer::cpu_timer;

const double p_kappa = 3.3;
const double p_lambda = 4.7;
const double p_beta = 1.0;

int main( int argc , char* argv[] )
{
    size_t N = 1024;
    size_t blocks = omp_get_max_threads();
    size_t steps = 100;
    size_t repeat = 5;
    bool split_range = true;
    if( argc > 1 ) N = boost::lexical_cast<size_t>( argv[1] );
    if( argc > 2 ) blocks = boost::lexical_cast<size_t>( argv[2] );
    if( argc > 3 ) steps = boost::lexical_cast<size_t>( argv[3] );
    if( argc > 4 ) repeat = boost::lexical_cast<size_t>( argv[4] );
    if( argc > 5 ) split_range = boost::lexical_cast<bool>( argv[5] );

    cout << "Size: " << N << " with " << blocks << " blocks and " << steps << " steps." << endl;

    accumulator_set< double, stats<tag::mean, tag::median> > acc_time;

    for(size_t n_rep = 0 ; n_rep != repeat ; n_rep++)
    {
        osc_chain system( p_kappa , p_lambda , p_beta );

        // fully random data
        vector<double> p_init( N ), q_init( N, 0 );
        uniform_real_distribution<double> distribution( 0.0 );
        mt19937 engine( 0 );
        auto generator = bind( distribution , engine );
        generate( p_init.begin() , p_init.end() , generator );

        if(split_range) {
            typedef openmp_state<double> state_type;
            typedef symplectic_rkn_sb3a_mclachlan<
                    state_type , state_type , double
                > stepper_type;

            omp_set_num_threads(blocks);
        
            // split into blocks
            state_type p( blocks );
            split(p_init, p);

            state_type q( blocks );
            split(q_init, q);

            clog << "split " << N << " into";
            for(size_t i = 0 ; i != p.size() ; i++)
                clog << ' ' << p[i].size();
            clog << endl;

            for(size_t n_run = 0 ; n_run != 5 ; n_run++) {
                cpu_timer timer;
                integrate_n_steps( stepper_type() , system ,
                                   make_pair( boost::ref(q) , boost::ref(p) ) ,
                                   0.0 , 0.01 , steps );
                double run_time = static_cast<double>(timer.elapsed().wall) * 1.0e-9;
                acc_time(run_time);
                clog << "run " << n_rep << "-" << n_run << " wall[s]: " << run_time << endl;
            }

        } else {
            typedef vector<double> state_type;
            typedef symplectic_rkn_sb3a_mclachlan<
                    state_type , state_type , double ,
                    state_type , state_type , double ,
                    openmp_range_algebra
                > stepper_type;
            
            omp_set_num_threads(blocks);

            state_type p(p_init), q(q_init);

            for(size_t n_run = 0 ; n_run != 5 ; n_run++) {
                cpu_timer timer;
                integrate_n_steps( stepper_type() , system ,
                                   make_pair( boost::ref(q) , boost::ref(p) ) ,
                                   0.0 , 0.01 , steps );
                double run_time = static_cast<double>(timer.elapsed().wall) * 1.0e-9;
                acc_time(run_time);
                clog << "run " << n_rep << "-" << n_run << " wall[s]: " << run_time << endl;
            }

        }
    }

    cout << " mean[s]: " << mean(acc_time)
         << " median[s]: " << median(acc_time) << endl;

    return 0;
}

#else

int main(int, char**) {
    return -1;
}

#endif
