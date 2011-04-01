/*
 * phase_osc_chain.cu
 *
 *  Created on: Apr 1, 2011
 *      Author: mario
 */

/*
 * This example shows how to use odeint on CUDA devices with thrust.
 * Note that we require at least Version 3.2 of the nVidia CUDA SDK
 * and the thrust library should be installed in the CUDA include
 * folder.
 *
 * As example we use a chain of phase oscillators with nearest neighbour
 * coupling, as described in:
 *
 * Avis H. Cohen, Philip J. Holmes and Richard H. Rand:
 * JOURNAL OF MATHEMATICAL BIOLOGY Volume 13, Number 3, 345-369,
 *
 */

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include <thrust/device_vector.h>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_operations.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_resize.hpp>

using namespace std;

using namespace boost::numeric::odeint;

//change this to float if your device does not support double computation
typedef double value_type;

//change this to host_vector< ... > of you want to run on CPU
typedef thrust::device_vector< value_type > state_type;
typedef thrust::device_vector< size_t > index_vector_type;
//typedef thrust::host_vector< value_type > state_type;
//typedef thrust::host_vector< size_t > index_vector_type;


/*
 * This implements the rhs of the dynamical equation:
 * \phi'_0 = \omega_0 + sin( \phi_1 - \phi_0 )
 * \phi'_i  = \omega_i + sin( \phi_i+1 - \phi_i ) + sin( \phi_i - \phi_i-1 )
 * \phi'_N-1 = \omega_N-1 + sin( \phi_N-1 - \phi_N-2 )
 */
class phase_oscillators
{

public:

    struct sys_functor
    {
        template< class Tuple >
        __host__ __device__
        void operator()( Tuple t )
        {
            const value_type phi = thrust::get<0>(t);
            const value_type phi_left = thrust::get<1>(t);
            const value_type phi_right = thrust::get<2>(t);
            const value_type omega = thrust::get<3>(t);
            // the dynamical equation
            thrust::get<4>(t) = omega + sin( phi_right - phi ) + sin( phi - phi_left );
        }
    };

    phase_oscillators( state_type &omega )
        : m_omega( omega ) , m_N( omega.size() ) , m_prev( m_N ) , m_next( m_N )
    {
        // build indices pointing to left and right neighbours
        thrust::counting_iterator<size_t> c( 0 );
        thrust::copy( c , c+m_N-1 , m_prev.begin()+1 );
        m_prev[0] = 0; // m_prev = { 0 , 0 , 1 , 2 , 3 , ... , N-1 }

        thrust::copy( c+1 , c+m_N , m_next.begin() );
        m_next[m_N-1] = m_N-1; // m_next = { 1 , 2 , 3 , ... , N-1 , N-1 }

        /*thrust::copy( m_prev.begin() , m_prev.end() ,
                    std::ostream_iterator< size_t >(std::cout, " ") );
        std::cout << std::endl;*/
    }



    void operator() ( const state_type &x , state_type &dxdt , const value_type dt )
    {
        thrust::for_each(
                thrust::make_zip_iterator(
                        thrust::make_tuple(
                                x.begin() ,
                                thrust::make_permutation_iterator( x.begin() , m_prev.begin() ) ,
                                thrust::make_permutation_iterator( x.begin() , m_next.begin() ) ,
                                m_omega.begin() ,
                                dxdt.begin()
                                ) ),
                thrust::make_zip_iterator(
                        thrust::make_tuple(
                                x.end() ,
                                thrust::make_permutation_iterator( x.begin() , m_prev.end() ) ,
                                thrust::make_permutation_iterator( x.begin() , m_next.end() ) ,
                                m_omega.end() ,
                                dxdt.end()) ) ,
                sys_functor()
                );
    }

private:
    const state_type &m_omega;
    const size_t m_N;
    index_vector_type m_prev;
    index_vector_type m_next;
};


const size_t N = 16;
const value_type epsilon = 6.0/(N*N); // should be < 8/N^2 to see phase locking

int main( int arc , char* argv[] )
{
    srand( time(NULL) );
    // create initial conditions on host:
    vector< value_type > x_host( N );
    //create omegas on host
    vector< value_type > omega_host( N );
    for( size_t i=0 ; i<N ; ++i )
    {
        x_host[i] = 2.0*3.14159265*(double)(rand())/RAND_MAX;
        omega_host[i] = (N-i)*epsilon; // decreasing frequencies
    }

    //copy to device
    state_type x = x_host;
    state_type omega = omega_host;

    //create error stepper
    explicit_rk4< state_type , value_type , state_type , value_type ,
                  thrust_algebra , thrust_operations , adjust_size_initially_tag  > stepper;

    phase_oscillators sys( omega );

    value_type t = 0.0;
    const value_type dt = 0.1;
    while( t < 10.0 )
    {
        stepper.do_step( sys , x , t , dt );
        t += dt;
    }

    /**ToDo: use integrate functions, maybe with algebra_dispatcher */

    //perform integration using standard Runge-Kutta-Cash-Carp Stepper and error bounds ~ 1E-6
    //integrate_const( phase_oscillators(omega) , x , 0.0 , 100.0 , 0.1 );

    thrust::copy( x.begin() , x.end() ,
            std::ostream_iterator< value_type >(std::cout, " ") );
    std::cout << std::endl;
}
