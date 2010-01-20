#include <iostream>

#include <boost/numeric/odeint.hpp>

//[ rhs_function
/* The type of container used to hold the state vector */
typedef std::vector<double> state_type;

const double gam = 0.15;

/* The rhs of x' = f(x) */
void harmonic_oscillator(const state_type &x, state_type &dxdt, const double t)
{
    dxdt[0] = x[1];
    dxdt[1] = -x[0] - gam*x[1];
}
//]

//[ rhs_class
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
//]

int main(int argc, char **argv)
{
    using namespace std;
    using namespace boost::numeric::odeint;

    //[ state_initialization
    state_type x(2);
    x[0] = 1.0; // start at x=1.0, p=0.0
    x[1] = 0.0;
    //]

    //[ integration
    vector<double> times;
    vector<state_type> x_t_vec;

    size_t steps = integrate( harmonic_oscillator , 
                              x , 0.0 , 10.0 , 
                              back_inserter( times ) ,
                              back_inserter( x_t_vec ) );
    //]

    /* the same as above using the class */
    /*
    harm_osc ho(0.15);
    steps = integrate( ho , 
                       x , 0.0 , 10.0 , 
                       back_inserter( times ) ,
                       back_inserter( x_t_vec ) );
    */

    //[ output
    for( size_t i=0; i<=steps; i++ ) { //initial state is 0th entry
        cout << times[i] << '\t' << x_t_vec[i][0] << '\t' << x_t_vec[i][1] << '\n';
    }
    //]
}
