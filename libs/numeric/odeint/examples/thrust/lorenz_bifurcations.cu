/*
 * phase_oscillator_ensemble.cu
 *
 * The example how the phase_oscillator ensemble can be implemented using CUDA and thrust
 *
 *  Created on: July 15, 2011
 *      Author: karsten
 */


#include <iostream>
#include <cmath>
#include <utility>

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_operations.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_resize.hpp>

#include <boost/random.hpp>



using namespace std;

using namespace boost::numeric::odeint;

//change this to float if your device does not support double computation
typedef double value_type;

//change this to host_vector< ... > of you want to run on CPU
typedef thrust::device_vector< value_type > state_type;
typedef thrust::device_vector< size_t > index_vector_type;
// typedef thrust::host_vector< value_type > state_type;
// typedef thrust::host_vector< size_t > index_vector_type;

const size_t N = 16384;
const value_type pi = 3.1415926535897932384626433832795029;


int main( int arc , char* argv[] )
{
    boost::mt19937 rng;
    boost::uniform_real< value_type > unif( 0.0 , 1.0 );
    boost::variate_generator< boost::mt19937&, boost::uniform_real< value_type > > gen( rng , unif );

    // vectors for host and device
    vector< value_type > x_host( N );
    state_type x( N );


    //create error stepper
    runge_kutta4< state_type , value_type , state_type , value_type , thrust_algebra , thrust_operations > stepper;
}
