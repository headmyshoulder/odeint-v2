/*
 * rk_test.cpp
 *
 *  Created on: Nov 23, 2010
 *      Author: mario
 */

#include <iostream>
#include <vector>
#include <tr1/array>
#include <boost/array.hpp>
#include <boost/fusion/container.hpp>

#include "runge_kutta_stepper.hpp"

using namespace std;

typedef tr1::array< double , 3 > state_type;
typedef runge_kutta_stepper< state_type , 1 > euler_stepper;
typedef runge_kutta_stepper< state_type , 2 > midpoint_stepper;


const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

const double dt = 0.001;

int main( void )
{

    boost::array< double , 1 > b = {{ 1.0 }};
    boost::array< double , 1 > c = {{ 0.0 }};

    euler_stepper euler( fusion::vector0<>() , b , c  );
    euler.print_vals();
    cout << typeid(euler_stepper::stage_vector_base).name() << endl;

//    boost::array< double , 1 > a2 = {{ 0.5 }};
//    boost::array< double , 2 > b2 = {{ 0.0 , 1.0 }};
//    boost::array< double , 2 > c2 = {{ 0.0 , 0.5 }};
//
//    cout << typeid(midpoint_stepper::coef_a_type).name() << endl;
//    cout << typeid(midpoint_stepper::coef_b_type).name() << endl;
//    cout << typeid(midpoint_stepper::coef_c_type).name() << endl;
//    cout << typeid(midpoint_stepper::stage_vector_base).name() << endl;
//
//    midpoint_stepper midpoint( fusion::make_vector( a2 ) , b2 , c2  );
//    midpoint.print_vals();

}
