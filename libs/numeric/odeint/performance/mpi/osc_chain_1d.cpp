/* Boost libs/numeric/odeint/performance/openmp/osc_chain_1d.cpp

 Copyright 2009-2012 Karsten Ahnert
 Copyright 2009-2012 Mario Mulansky

 stronlgy nonlinear hamiltonian lattice in 2d

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <iostream>
#include <vector>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/mpi/mpi.hpp>

#include <boost/program_options.hpp>
#include <boost/random.hpp>
#include <boost/timer/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include "osc_chain_1d_system.hpp"

using namespace std;
using namespace boost::numeric::odeint;
using namespace boost::accumulators;
using namespace boost::random;
using namespace boost::program_options;

using boost::timer::cpu_timer;

const double p_kappa = 3.3;
const double p_lambda = 4.7;
const double p_beta = 1.0;

typedef vector<double> inner_state_type;
typedef mpi_state< inner_state_type > state_type;
typedef symplectic_rkn_sb3a_mclachlan<
        state_type , state_type , double
    > stepper_type;


int main( int argc , char* argv[] )
{
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    size_t N, steps, repeat;
    bool dump;
    options_description desc("Options");
    desc.add_options()
        ("help,h", "show this help")
        ("length", value(&N)->default_value(1024), "length of chain")
        ("steps", value(&steps)->default_value(100), "simulation steps")
        ("repeat", value(&repeat)->default_value(25), "repeat runs")
        ("dump", bool_switch(&dump), "dump final state to stderr (on node 0)")
        ;
    variables_map vm;
    store(command_line_parser(argc, argv).options(desc).run(), vm);
    notify(vm);
    if(vm.count("help"))
    {
        if(world.rank() == 0)
            cerr << desc << endl;
        return EXIT_FAILURE;
    }
    cout << "length\tsteps\tthreads\ttime" << endl;

    accumulator_set< double, stats<tag::mean, tag::median> > acc_time;

    for(size_t n_rep = 0 ; n_rep != repeat ; n_rep++)
    {
        osc_chain system( p_kappa , p_lambda , p_beta );

        // fully random data
        inner_state_type p_init( N ), q_init( N, 0 );
        if(world.rank() == 0) {
            uniform_real_distribution<double> distribution( 0.0 );
            mt19937 engine( 0 );
            generate( p_init.begin() , p_init.end() , boost::bind( distribution , engine ) );
        }

        // send to nodes
        state_type p( world ), q( world );
        split(p_init, p);
        split(q_init, q);

        for(size_t n_run = 0 ; n_run != 5 ; n_run++) {
            cpu_timer timer;
            world.barrier();
            integrate_n_steps( stepper_type() , system ,
                               make_pair( boost::ref(q) , boost::ref(p) ) ,
                               0.0 , 0.01 , steps );
            world.barrier();
            if(world.rank() == 0) {
                double run_time = static_cast<double>(timer.elapsed().wall) * 1.0e-9;
                acc_time(run_time);
                clog << "run " << n_rep << "-" << n_run << " wall[s]: " << run_time << endl;
            }
        }
    }

    if(world.rank() == 0)
        cout << " mean[s]: " << mean(acc_time)
             << " median[s]: " << median(acc_time) << endl;

    return 0;
}

