/* Boost numeric/odeint/examples/lorenz_stepper.cpp
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 A simple compile test program

 Furthermore, the usage of std::tr1::array and std::vector in odeint is
 shown and the performance of both containers is compared.

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <iostream>
#include <vector>
#include <list>
#include <tr1/array>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/container_traits_blitz_array.hpp>
#include <boost/numeric/odeint/container_traits_mtl4_dense2d.hpp>
#include <boost/numeric/odeint/container_traits_ublas_matrix.hpp>

#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;

const size_t n1 = 3;
const size_t n2 = 5;
const double dt = 0.01;

typedef std::vector< double > state_type1;
typedef std::tr1::array< double , n1 > state_type2;
typedef blitz::Array< double , 2 > state_type3;
typedef mtl::dense2D< double > state_type4;
typedef boost::numeric::ublas::matrix< double > state_type5;

void deriv1( state_type1 &x , state_type1 &dxdt , double t )
{
    const double sigma = 10.0 , R = 28.0 , b = 8.0 / 3.0;
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

void deriv2( state_type2 &x , state_type2 &dxdt , double t )
{
    const double sigma = 10.0 , R = 28.0 , b = 8.0 / 3.0;
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

void deriv3( state_type3 &x , state_type3 &dxdt , double t )
{
    for( int i=0 ; i<int(n1) ; ++i )
        for( int j=0 ; j<int(n2) ; ++j ) 
            dxdt( i , j ) = drand48();
}

void deriv4( state_type4 &x , state_type4 &dxdt , double t )
{
    for( size_t i=0 ; i<n1 ; ++i )
        for( size_t j=0 ; j<n2 ; ++j ) 
            dxdt( i , j ) = drand48();
}

void deriv5( state_type5 &x , state_type5 &dxdt , double t )
{
    for( size_t i=0 ; i<n1 ; ++i )
        for( size_t j=0 ; j<n2 ; ++j ) 
            dxdt( i , j ) = drand48();
}

template<
    class Stepper1 ,
    class Stepper2 ,
    class Stepper3 ,
    class Stepper4 ,
    class Stepper5
    >
void test_steppers(
    Stepper1 &stepper1 ,
    Stepper2 &stepper2 ,
    Stepper3 &stepper3 ,
    Stepper4 &stepper4 ,
    Stepper5 &stepper5
    )
{
    state_type1 x1( n1 , 0.0 ) , dxdt1( x1 );
    state_type2 x2 = {{ 1.0 , 0.0 , 0.0 }} , dxdt2( x2 );
    state_type3 x3( n1 , n2 ) , dxdt3( x3 );
    state_type4 x4( n1 , n2 ) , dxdt4( x4 );
    state_type5 x5( n1 , n2 ) , dxdt5( x5 );

    stepper1.do_step( deriv1 , x1 , 0.0 , dt );
    stepper2.do_step( deriv2 , x2 , 0.0 , dt );
    stepper3.do_step( deriv3 , x3 , 0.0 , dt );
    stepper4.do_step( deriv4 , x4 , 0.0 , dt );
    stepper5.do_step( deriv5 , x5 , 0.0 , dt );

    stepper1.do_step( deriv1 , x1 , dxdt1 , 0.0 , dt );
    stepper2.do_step( deriv2 , x2 , dxdt2 , 0.0 , dt );
    stepper3.do_step( deriv3 , x3 , dxdt3 , 0.0 , dt );
    stepper4.do_step( deriv4 , x4 , dxdt4 , 0.0 , dt );
    stepper5.do_step( deriv5 , x5 , dxdt5 , 0.0 , dt );
}



template<
    class Stepper1 ,
    class Stepper2 ,
    class Stepper3 ,
    class Stepper4 ,
    class Stepper5
    >
void test_error_steppers(
    Stepper1 &stepper1 ,
    Stepper2 &stepper2 ,
    Stepper3 &stepper3 ,
    Stepper4 &stepper4 ,
    Stepper5 &stepper5
    )
{
    state_type1 x1( n1 , 0.0 ) , dxdt1( x1 ) , error1( x1 );
    state_type2 x2 = {{ 1.0 , 0.0 , 0.0 }} , dxdt2( x2 ) , error2( x2 );
    state_type3 x3( n1 , n2 ) , dxdt3( x3 ) , error3( x3 );
    state_type4 x4( n1 , n2 ) , dxdt4( x4 ) , error4( x4 );
    state_type5 x5( n1 , n2 ) , dxdt5( x5 ) , error5( x5 );

    stepper1.do_step( deriv1 , x1 , 0.0 , dt , error1 );
    stepper2.do_step( deriv2 , x2 , 0.0 , dt , error2 );
    stepper3.do_step( deriv3 , x3 , 0.0 , dt , error3 );
    stepper4.do_step( deriv4 , x4 , 0.0 , dt , error4 );
    stepper5.do_step( deriv5 , x5 , 0.0 , dt , error5 );

    stepper1.do_step( deriv1 , x1 , dxdt1 , 0.0 , dt , error1 );
    stepper2.do_step( deriv2 , x2 , dxdt2 , 0.0 , dt , error2 );
    stepper3.do_step( deriv3 , x3 , dxdt3 , 0.0 , dt , error3 );
    stepper4.do_step( deriv4 , x4 , dxdt4 , 0.0 , dt , error4 );
    stepper5.do_step( deriv5 , x5 , dxdt5 , 0.0 , dt , error5 );
}

template< class Stepper , class Deriv , class Container >
void test_stepper( Stepper &stepper , Deriv &deriv , Container &x , Container &dxdt )
{
    typename Stepper::container_type xx;
    typename Stepper::order_type order;
    typename Stepper::time_type time;
    typename Stepper::traits_type traits;
    typename Stepper::value_type value;
    typename Stepper::iterator iterator;
    typename Stepper::const_iterator const_iteratio;

    stepper.do_step( deriv , x , 0.0 , dt );
    stepper.do_step( deriv , x , dxdt , 0.0 , dt );
}



int main( int argc , char **argv )
{
    cout << "Hello world!" << endl;

    state_type1 x1( n1 , 0.0 ) , dxdt1( x1 ) , error1( x1 );
    state_type2 x2 = {{ 1.0 , 0.0 , 0.0 }} , dxdt2( x2 ) , error2( x2 );
    state_type3 x3( n1 , n2 ) , dxdt3( x3 ) , error3( x3 );
    state_type4 x4( n1 , n2 ) , dxdt4( x4 ) , error4( x4 );
    state_type5 x5( n1 , n2 ) , dxdt5( x5 ) , error5( x5 );

    stepper_euler< state_type1 > s1;
    stepper_euler< state_type2 > s2;
    stepper_euler< state_type3 > s3;
    stepper_euler< state_type4 > s4;
    stepper_euler< state_type5 > s5;

    test_stepper( s1 , deriv1 , x1 , dxdt1 );
    

    {
        stepper_euler< state_type1 > s1;
        stepper_euler< state_type2 > s2;
        stepper_euler< state_type3 > s3;
        stepper_euler< state_type4 > s4;
        stepper_euler< state_type5 > s5;
        test_steppers( s1 , s2 , s3 , s4 , s5 );
    }
    cout << "Euler Stepper OK!" << endl;


    {
        stepper_half_step< stepper_euler< state_type1 > > s1;
        stepper_half_step< stepper_euler< state_type2 > > s2;
        stepper_half_step< stepper_euler< state_type3 > > s3;
        stepper_half_step< stepper_euler< state_type4 > > s4;
        stepper_half_step< stepper_euler< state_type5 > > s5;
        test_steppers( s1 , s2 , s3 , s4 , s5 );
        test_error_steppers( s1 , s2 , s3 , s4 , s5 );
    }
    cout << "Half Step Stepper OK!" << endl;


    {
        stepper_rk4< state_type1 > s1;
        stepper_rk4< state_type2 > s2;
        stepper_rk4< state_type3 > s3;
        stepper_rk4< state_type4 > s4;
        stepper_rk4< state_type5 > s5;
        test_steppers( s1 , s2 , s3 , s4 , s5 );
    }
    cout << "RK4 Stepper OK!" << endl;




    {
        stepper_rk4_classical< state_type1 > s1;
        stepper_rk4_classical< state_type2 > s2;
        stepper_rk4_classical< state_type3 > s3;
        stepper_rk4_classical< state_type4 > s4;
        stepper_rk4_classical< state_type5 > s5;
        test_steppers( s1 , s2 , s3 , s4 , s5 );
    }
    cout << "RK4 Classical Stepper OK!" << endl;



    {
        stepper_rk5_ck< state_type1 > s1;
        stepper_rk5_ck< state_type2 > s2;
        stepper_rk5_ck< state_type3 > s3;
        stepper_rk5_ck< state_type4 > s4;
        stepper_rk5_ck< state_type5 > s5;
        test_error_steppers( s1 , s2 , s3 , s4 , s5 );
    }
    cout << "RK5 CK Stepper OK!" << endl;


    {
        stepper_midpoint< state_type1 > s1;
        stepper_midpoint< state_type2 > s2;
        stepper_midpoint< state_type3 > s3;
        stepper_midpoint< state_type4 > s4;
        stepper_midpoint< state_type5 > s5;
        test_steppers( s1 , s2 , s3 , s4 , s5 );
    }
    cout << "Midpoint Stepper OK!" << endl;


    {
        stepper_rk78_fehlberg< state_type1 > s1;
        stepper_rk78_fehlberg< state_type2 > s2;
        stepper_rk78_fehlberg< state_type3 > s3;
        stepper_rk78_fehlberg< state_type4 > s4;
        stepper_rk78_fehlberg< state_type5 > s5;
        test_steppers( s1 , s2 , s3 , s4 , s5 );
    }
    cout << "RK78 Stepper OK!" << endl;


    cout << "Everything is allright!" << endl;

    return 0;
}
