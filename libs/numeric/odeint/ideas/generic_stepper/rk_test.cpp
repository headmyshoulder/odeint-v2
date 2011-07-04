#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint/stepper/explicit_rk4.hpp>

#include "runge_kutta_stepper.hpp"

using namespace std;

namespace odeint = boost::numeric::odeint;

typedef boost::array< double , 3 > state_type;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

int main( int argc , char **argv )
{
	typedef runge_kutta_stepper< state_type , 4 > rk_type;
	typedef rk_type::coef_a_type coef_a_type;
	typedef rk_type::coef_b_type coef_b_type;
	typedef rk_type::coef_c_type coef_c_type;

	const boost::array< double , 1 > a1 = {{ 0.5 }};
	const boost::array< double , 2 > a2 = {{ 0.0 , 0.5 }};
	const boost::array< double , 3 > a3 = {{ 0.0 , 0.0 , 1.0 }};

	const coef_a_type a = fusion::make_vector( a1 , a2 , a3 );
	const coef_b_type b = {{ 1.0/6 , 1.0/3 , 1.0/3 , 1.0/6 }};
	const coef_c_type c = {{ 0.0 , 0.5 , 0.5 , 1.0 }};

	rk_type rk( a , b , c );
	rk.print_vals();

	state_type x = {{ 1.0 , 1.0 , 2.0 }};
	double t = 0.0;
	rk.do_step( lorenz , x , t , 0.1 );

	cout.precision(16);

	cout << x[0] << " , " << x[1] << " , " << x[2] << endl;


	odeint::explicit_rk4< state_type > rk4;
	state_type x2 = {{ 1.0 , 1.0 , 2.0 }};
	rk4.do_step( lorenz , x2 , t , 0.1 );
	cout << x2[0] << " , " << x2[1] << " , " << x2[2] << endl;
}
