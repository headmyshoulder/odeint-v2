#include <iostream>
#include <cmath>

#include <boost/timer.hpp>

#include <boost/numeric/odeint/stepper/rosenbrock4_controller.hpp>
#include <boost/numeric/odeint/stepper/generic_dense_output_stepper.hpp>
#include <boost/numeric/odeint/stepper/generation.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>

#define tab "\t"

typedef boost::numeric::odeint::rosenbrock4< double > stepper_type;
typedef stepper_type::state_type state_type;
typedef stepper_type::matrix_type matrix_type;

const double mu = 1000.0;


struct vdp_stiff
{
    template< class State , class Deriv >
    void operator()( const State &x , Deriv &dxdt , double t )
    {
        dxdt[0] = x[1];
        dxdt[1] = -x[0] - mu * x[1] * (x[0]*x[0]-1.0);
    }
};

struct vdp_stiff_jacobi
{
    void operator()( const state_type &x , matrix_type &J , const double &t , state_type &dfdt )
    {
        J(0, 0) = 0.0;
        J(0, 1) = 1.0;
        J(1, 0) = -1.0 - 2.0*mu * x[0] * x[1];
        J(1, 1) = -mu * ( x[0] * x[0] - 1.0);

        dfdt[0] = 0.0;
        dfdt[1] = 0.0;
    }
};


int main( int argc , char **argv )
{
    using namespace std;
    using namespace boost::numeric::odeint;

    stepper_type stepper1;

    state_type x1( 2 );
    x1[ 0 ] = 1.0 ; x1[ 1 ] = 1.0;
    state_type x2 = x1;
    state_type x3 = x1;
    state_type x4 = x1;

    boost::timer timer;
    const double t_max = 10000000.0;

    timer.restart();
    size_t steps1 = integrate_adaptive( make_dense_output( 1.0e-6 , 1.0e-6 , stepper1 ) ,
            make_pair( vdp_stiff() , vdp_stiff_jacobi() ) , x1 , 0.0 , t_max , 0.1 );
    double t1 = timer.elapsed();
    cout << steps1 << tab << t1 << tab << x1[0] << tab << x1[1] << endl;

    generic_dense_output_stepper< rosenbrock4_controller< rosenbrock4< double > > > stepper2;
    timer.restart();
    size_t steps2 = integrate_adaptive( stepper2 , make_pair( vdp_stiff() , vdp_stiff_jacobi() ) , x2 , 0.0 , t_max , 0.1 );
    double t2 = timer.elapsed();
    cout << steps2 << tab << t2 << tab << x2[0] << tab << x2[1] << endl;



    return 0;
}

