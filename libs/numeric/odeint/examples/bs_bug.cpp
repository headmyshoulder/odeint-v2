#include <iostream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array< double , 3 > state_type;

typedef boost::numeric::odeint::bulirsch_stoer< state_type > bulirsch_stoer_type;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

int main() {
    
    const double eps_abs = 10.0;//1E-4;
    const double eps_rel = 1E-4;
    const double t_end = 10.0;

    bulirsch_stoer_type bulirsch_stoer( eps_abs , eps_rel );
    state_type x = {{ 10.0 , 10.0 , 10.0 }};

    cout.precision(17);

    size_t steps = integrate_adaptive( bulirsch_stoer , lorenz , x , 0.0 , t_end , 0.1 );

    cout << steps << " steps" << endl;

    return 0;
}
