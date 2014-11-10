/*
 * Simulation of an ensemble of Roessler attractors
 *
 * Copyright 2014 Mario Mulansky
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */


#include <iostream>
#include <vector>
#include <random>

#include <boost/timer.hpp>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

namespace odeint = boost::numeric::odeint;

typedef boost::timer timer_type;

typedef boost::array<double, 3> state_type;
typedef std::vector<state_type> state_vec;

//---------------------------------------------------------------------------
struct roessler_system {
    const double m_a, m_b, m_c;

    roessler_system(const double a, const double b, const double c)
        : m_a(a), m_b(b), m_c(c)
    {}

    void operator()(const state_type &x, state_type &dxdt, const double t) const
    {
        dxdt[0] = -x[1] - x[2];
        dxdt[1] = x[0] + m_a * x[1];
        dxdt[2] = m_b + x[2] * (x[0] - m_c);
    }
};

//---------------------------------------------------------------------------
int main(int argc, char *argv[]) {
if(argc<3)
{
    std::cerr << "Expected size and steps as parameter" << std::endl;
    exit(1);
}
const size_t n = atoi(argv[1]);
const size_t steps = atoi(argv[2]);
//const size_t steps = 50;

const double dt = 0.01;

const double a = 0.2;
const double b = 1.0;
const double c = 9.0;

// random initial conditions on the device
std::vector<double> x(n), y(n), z(n);
std::default_random_engine generator;
std::uniform_real_distribution<double> distribution_xy(-8.0, 8.0);
std::uniform_real_distribution<double> distribution_z(0.0, 20.0);
auto rand_xy = std::bind(distribution_xy, std::ref(generator));
auto rand_z = std::bind(distribution_z, std::ref(generator));
std::generate(x.begin(), x.end(), rand_xy);
std::generate(y.begin(), y.end(), rand_xy);
std::generate(z.begin(), z.end(), rand_z);

state_vec state(n);
for(size_t i=0; i<n; ++i)
{
    state[i][0] = x[i];
    state[i][1] = y[i];
    state[i][2] = z[i];
}

std::cout << "# n: " << n << std::endl;

// Stepper type - use never_resizer for slight performance improvement
odeint::runge_kutta4_classic<state_type, double, state_type, double,
                             odeint::array_algebra,
                             odeint::default_operations,
                             odeint::never_resizer> stepper;

roessler_system sys(a, b, c);

timer_type timer;

double t = 0.0;

for (int step = 0; step < steps; step++)
{
    for(size_t i=0; i<n; ++i)
    {
        stepper.do_step(sys, state[i], t, dt);
    }
    t += dt;
}

std::cout << "Integration finished, runtime for " << steps << " steps: ";
std::cout << timer.elapsed() << " s" << std::endl;

// compute some accumulation to make sure all results have been computed
double s = 0.0;
for(size_t i = 0; i < n; ++i)
{
    s += state[i][0];
}

std::cout << x[0] << std::endl;
std::cout << s/n << std::endl;

}
