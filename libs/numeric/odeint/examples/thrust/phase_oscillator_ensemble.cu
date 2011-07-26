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


struct calc_mean_field
{
    struct sin_functor : public thrust::unary_function< value_type , value_type >
    {
        __host__ __device__
        value_type operator()( value_type x) const
        {
            return sin( x );
        }
    };

    struct cos_functor : public thrust::unary_function< value_type , value_type >
    {
        __host__ __device__
        value_type operator()( value_type x) const
        {
            return cos( x );
        }
    };

    std::pair< value_type , value_type > get_mean( const state_type &x ) const
    {
        value_type sin_sum = thrust::reduce(
                thrust::make_transform_iterator( x.begin() , sin_functor() ) ,
                thrust::make_transform_iterator( x.end() , sin_functor() ) );
        value_type cos_sum = thrust::reduce(
                thrust::make_transform_iterator( x.begin() , cos_functor() ) ,
                thrust::make_transform_iterator( x.end() , cos_functor() ) );

        cos_sum /= value_type( x.size() );
        sin_sum /= value_type( x.size() );

        value_type K = sqrt( cos_sum * cos_sum + sin_sum * sin_sum );
        value_type Theta = atan2( sin_sum , cos_sum );

        return std::make_pair( K , Theta );
    }
};

class phase_oscillator_ensemble
{

public:

    struct sys_functor
    {
        value_type m_K , m_Theta , m_epsilon;
        sys_functor( value_type K , value_type Theta , value_type epsilon )
        : m_K( K ) , m_Theta( Theta ) , m_epsilon( epsilon ) { }

        template< class Tuple >
        __host__ __device__
        void operator()( Tuple t )
        {
            thrust::get<2>(t) = thrust::get<1>(t) + m_epsilon * m_K * sin( m_Theta - thrust::get<0>(t) );
        }
    };

    phase_oscillator_ensemble( size_t N , value_type g = 1.0 , value_type epsilon = 1.0 )
        : m_omega() , m_N( N ) , m_epsilon( epsilon )
    {
        create_frequencies( g );
    }

    void create_frequencies( value_type g )
    {
        boost::mt19937 rng;
        boost::cauchy_distribution< value_type > cauchy( 0.0 , g );
        boost::variate_generator< boost::mt19937&, boost::cauchy_distribution< value_type > > gen( rng , cauchy );
        vector< value_type > omega( m_N );
        generate( omega.begin() , omega.end() , gen );
        m_omega = omega;
    }

    void set_epsilon( value_type epsilon ) { m_epsilon = epsilon; }

    value_type get_epsilon( void ) const { return m_epsilon; }

    void operator() ( const state_type &x , state_type &dxdt , const value_type dt ) const
    {
        calc_mean_field mean_field_calculator;
        std::pair< value_type , value_type > mean_field = mean_field_calculator.get_mean( x );

        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple( x.begin() , m_omega.begin() , dxdt.begin() ) ),
                thrust::make_zip_iterator( thrust::make_tuple( x.end() , m_omega.end() , dxdt.end()) ) ,
                sys_functor( mean_field.first , mean_field.second , m_epsilon )
                );
    }

private:

    state_type m_omega;
    const size_t m_N;
    value_type m_epsilon;
};


//[ phase_oscillator_ensemble_observer
struct statistics_observer
{
    value_type m_K_mean;
    size_t m_count;

    statistics_observer( void )
    : m_K_mean( 0.0 ) , m_count( 0 ) { }

    template< class State >
    void operator()( const State &x , value_type t )
    {
        calc_mean_field mean_field_calculator;
        std::pair< value_type , value_type > mean = mean_field_calculator.get_mean( x );
        m_K_mean += mean.first;
        ++m_count;
    }

    value_type get_K_mean( void ) const { return ( m_count != 0 ) ? m_K_mean / value_type( m_count ) : 0.0 ; }

    void reset( void ) { m_K_mean = 0.0; m_count = 0; }
};
//]



// const size_t N = 16384 * 128;
const size_t N = 16384;
const value_type pi = 3.1415926535897932384626433832795029;
const value_type dt = 0.1;

int main( int arc , char* argv[] )
{
    boost::mt19937 rng;
    boost::uniform_real< value_type > unif( 0.0 , 2.0 * pi );
    boost::variate_generator< boost::mt19937&, boost::uniform_real< value_type > > gen( rng , unif );

    // vectors for host and device
    vector< value_type > x_host( N );
    state_type x( N );


    //create error stepper
    runge_kutta4< state_type , value_type , state_type , value_type , thrust_algebra , thrust_operations > stepper;
    phase_oscillator_ensemble ensemble( N , 1.0 );
    statistics_observer obs;

    for( value_type epsilon = 0.0 ; epsilon < 5.0 ; epsilon += 0.1 )
    {
        ensemble.set_epsilon( epsilon );
        obs.reset();

        // start with random initial conditions
        generate( x_host.begin() , x_host.end() , gen );
        x = x_host;

        // calculate some transients steps
        integrate_const( stepper , boost::ref( ensemble ) , x , 0.0 , 10.0 , dt );

        // integrate and compute the statistics
        integrate_const( stepper , boost::ref( ensemble ) , x , 0.0 , 100.0 , dt , boost::ref( obs ) );
        cout << epsilon << "\t" << obs.get_K_mean() << endl;
    }
}
