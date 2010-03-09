#include <iostream>

#include <boost/numeric/odeint.hpp>

using namespace std;

/* The type of container used to hold the state vector */
typedef std::vector<double> state_type;

class harm_osc {

    double m_gam;

public:
    harm_osc( double gam ) : m_gam(gam) { }

    void operator() (const state_type &x, state_type &dxdt, const double t)
    {
        dxdt[0] = x[1];
        dxdt[1] = -x[0] - m_gam*x[1];
    }
};

int main(int argc, char **argv)
{
    using namespace std;
    using namespace boost::numeric::odeint;

    state_type x(2);
    x[0] = 1.0; // start at x=1.0, p=0.0
    x[1] = 0.0;

    harm_osc harmonic_oscillator(0.15);

    //[ define_const_stepper
    stepper_rk4< state_type > rk4;
    //]

    //[ integrate_const
    integrate_const( rk4, harmonic_oscillator, x , 0.0, 10.0 , 0.01);
    //]


    //[ define_adapt_stepper
    stepper_rk5_ck< state_type > rk5;
    //]

    //[ define_conntrolled_stepper
    controlled_stepper_standard< stepper_rk5_ck< state_type > > 
        controlled_rk5( 1E-6 , 1E-7 , 1.0 , 1.0 );
    //]

    //[ integrate_adapt
    integrate_adaptive( controlled_rk5 ,
                        harmonic_oscillator, x, 0.0, 10.0, 0.01  );
    //]
    
    cout << x[0] << '\t' << x[1] << endl;

}
