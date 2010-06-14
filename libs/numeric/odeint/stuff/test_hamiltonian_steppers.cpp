#include <iostream>
#include <vector>
#include <utility>
#include <boost/numeric/odeint.hpp>

#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;

typedef vector< double > coord_type;
typedef pair< coord_type , coord_type > state_type;
typedef hamiltonian_stepper_euler< coord_type > euler_stepper_type;
typedef hamiltonian_stepper_euler_qfunc< coord_type > euler_stepper_type_qfunc;
typedef hamiltonian_stepper_rk< coord_type > rk_stepper_type;
typedef hamiltonian_stepper_rk_qfunc< coord_type > rk_stepper_type_qfunc;

static const double omega0 = 1.0;
static const double omega1 = 2.0;

void qfunc(const coord_type &q, coord_type &dpdt )
{
    dpdt[0] = -omega0*omega0*q[0];// + q[1];
    dpdt[1] = -omega1*omega1*q[1];// + q[0];
}

void pfunc(const coord_type &p, coord_type &dqdt )
{
    dqdt[0] = p[0];
    dqdt[1] = p[1];
}

double energy( state_type &state )
{
    coord_type &q = state.first;
    coord_type &p = state.second;
    return p[0]*p[0]/2.0 + p[1]*p[1]/2.0 
        + omega0*omega0*q[0]*q[0]/2.0 
        + omega1*omega1*q[1]*q[1]/2.0;
    //- q[0]*q[1];
}


int main( int argc , char **argv )
{
    coord_type x1(2) , p1(2 , 0.0);
    x1[0] = 1.0;
    x1[1] = -1.0;
    coord_type x2(x1) , x3(x1) , x4(x1);
    coord_type p2(p1) , p3(p1) , p4(p1);
    state_type state1 = make_pair( x1 , p1 );
    state_type state2 = make_pair( x2 , p2 );
    state_type state3 = make_pair( x3 , p3 );
    state_type state4 = make_pair( x4 , p4 );
    double initial_energy = energy( state1 );
    clog << "initial energy: " << initial_energy << endl;

    euler_stepper_type euler;
    euler_stepper_type_qfunc euler_qfunc;
    rk_stepper_type rk;
    rk_stepper_type_qfunc rk_qfunc;

    const double dt = 0.1;
    double t = 0.0;
    cout.precision(16);
    for( int i=0 ; i<100 ; ++i , t += dt)
    {
        euler.do_step( make_pair( &pfunc , &qfunc ) , state1 , t , dt );
        euler_qfunc.do_step( &qfunc , state2 , t , dt );

        rk.do_step( make_pair( &pfunc , &qfunc ) , state3 , t , dt );
        rk_qfunc.do_step( &qfunc , state4 , t , dt );

        cout << t << tab << energy( state1 ) - initial_energy << tab;
        cout << energy( state2 ) - initial_energy << tab;
        cout << energy( state3 ) - initial_energy << tab;
        cout << energy( state4 ) - initial_energy << endl;
    }
    clog << "final energy : " << energy(state1) << tab;
    clog << energy( state2 ) << tab;
    clog << energy( state3 ) << tab;
    clog << energy( state4 ) << endl;

}
