#include <iostream>
#include <vector>

#include <boost/numeric/odeint.hpp>



//[ rhs_function
/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

const double gam = 0.15;

/* The rhs of x' = f(x) */
void harmonic_oscillator( const state_type &x , state_type &dxdt , const double t )
{
    dxdt[0] = x[1];
    dxdt[1] = -x[0] - gam*x[1];
}
//]





//[ rhs_class
/* The rhs of x' = f(x) defined as a class */
class harm_osc {

    double m_gam;

public:
    harm_osc( double gam ) : m_gam(gam) { }

    void operator() ( const state_type &x , state_type &dxdt , const double t )
    {
        dxdt[0] = x[1];
        dxdt[1] = -x[0] - m_gam*x[1];
    }
};
//]





//[ integrate_observer
struct push_back_state_and_time
{
	std::vector< state_type >& m_states;
	std::vector< double >& m_times;

	push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
	: m_states( states ) , m_times( times ) { }

	void operator()( const state_type &x , double t )
	{
		m_states.push_back( x );
		m_times.push_back( t );
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
    size_t steps = integrate( harmonic_oscillator ,
                              x , 0.0 , 10.0 , 0.1 );
    //]



    //[ integration_class
    harm_osc ho(0.15);
    steps = integrate( ho ,
                       x , 0.0 , 10.0 , 0.1 );
    //]





    //[ integrate_observ
    vector<state_type> x_vec;
    vector<double> times;

    steps = integrate( harmonic_oscillator ,
                       x , 0.0 , 10.0 , 0.1 ,
                       push_back_state_and_time( x_vec , times ) );

    /* output */
    for( size_t i=0; i<=steps; i++ )
    {
        cout << times[i] << '\t' << x_vec[i][0] << '\t' << x_vec[i][1] << '\n';
    }
    //]







    //[ define_const_stepper
    explicit_rk4< state_type > stepper;
    integrate_const( stepper , harmonic_oscillator , x , 0.0 , 10.0 , 0.01 );
    //]




    //[ integrate_const_loop
    const double dt = 0.01;
    for( double t=0.0 ; t<10.0 ; t+= dt )
    	stepper.do_step( harmonic_oscillator , x , t , dt );
    //]




    //[ define_adapt_stepper
    typedef explicit_error_rk54_ck< state_type > error_stepper_type;
    error_stepper_type rk54;
    //]



    //[ integrate_adapt
    typedef controlled_error_stepper< error_stepper_type > controlled_stepper_type;
    controlled_stepper_type controlled_stepper;
    integrate_adaptive( controlled_stepper , harmonic_oscillator , x , 0.0 , 10.0 , 0.01 );
    //]
}
