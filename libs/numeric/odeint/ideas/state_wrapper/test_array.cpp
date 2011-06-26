#include <iostream>
#include <boost/array.hpp>

#include "explicit_euler.hpp"
#include "size_adjuster.hpp"

using namespace std;

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

int main()
{
    explicit_euler< state_type , initially_resizer > euler;
    //you would use never_resizer for arrays, but just to show that no resizing is called
    //even with initially_resizer I use this one here.

    state_type x;
    x[0] = 1.0; x[1] = 1.0; x[2] = 2.0;

    euler.do_step( lorenz , x , 0.0 , 0.1 );

    cout << x[0] << "  " << x[1] << "  " << x[2] << endl;
}
