#include <iostream>
#include <vector>

#define DECORATE_CALLS
#include <boost/numeric/odeint.hpp>
#include <thrust/device_vector.h>
#include <thrust/for_each.h>
#include <thrust/iterator/zip_iterator.h>

#include <stdio.h>
#include <cuda_runtime.h>

namespace odeint = boost::numeric::odeint;



//---------------------------------------------------------------------------
template <typename T>
struct point3d {
    T x;
    T y;
    T z;
};

template <typename T>
__host__ __device__ point3d<T> operator+(point3d<T> a, point3d<T> b) {
    point3d<T> c = {a.x + b.x, a.y + b.y, a.z + b.z};
    return c;
}

template< typename T>
__host__ __device__ point3d<T> operator+(T a, point3d<T> b)
{
    point3d<T> c={a + b.x, a + b.y, a + b.z};
    return c;
}

template< typename T>
__host__ __device__ point3d<T> operator+(point3d<T> b, T a)
{
    point3d<T> c={a + b.x, a + b.y, a + b.z};
    return c;
}


template <typename T>
__host__ __device__ point3d<T> operator-(point3d<T> a, point3d<T> b) {
    point3d<T> c = {a.x - b.x, a.y - b.y, a.z - b.z};
    return c;
}

template <typename T>
__host__ __device__ point3d<T> operator*(T a, point3d<T> b) {
    point3d<T> c = {a * b.x, a * b.y, a * b.z};
    return c;
}

template <typename T>
__host__ __device__ point3d<T> operator/(point3d<T> a, point3d<T> b) {
    point3d<T> c = {a.x / b.x, a.y / b.y, a.z / b.z};
    return c;
}


template <typename T>
__host__ __device__ point3d<T> abs(point3d<T> p)
{
    point3d<T> ret;
    ret.x = abs( p.x );
    ret.y = abs( p.y );
    ret.z = abs( p.z );
    return ret;
}

namespace boost { namespace numeric { namespace odeint { 
template<typename T>
struct is_resizeable< point3d<T> > : boost::false_type { };
} } } 

namespace boost { namespace numeric { namespace odeint { 
template<typename T>
struct vector_space_reduce< point3d< T > >
{
    template< class Op >
    __host__ __device__ T operator()( const point3d<T> &x , Op op , T init ) const
    {
        init = op( init , x.x );
        init = op( init , x.y );
        init = op( init , x.z );
        return init;
    }
};
} } } 
//---------------------------------------------------------------------------



typedef double value_type;
typedef point3d<value_type> state_type;

//---------------------------------------------------------------------------
std::ostream& operator<<(std::ostream &os, const state_type &s) {
    return os << "[" << s.x << " " << s.y << " " << s.z << "]";
}


const value_type sigma = 10.0;
const value_type b = 8.0 / 3.0;
const value_type t_max = 10.0;

struct lorenz_system {
    value_type R;
    
    lorenz_system(value_type r = 0) : R(r) {}

    __host__ __device__ void operator()(const state_type &s, state_type &dsdt, value_type t) {
	dsdt.x = sigma * (s.y - s.x);
	dsdt.y = R * s.x - s.y - s.x * s.z;
	dsdt.z = s.x * s.y - b * s.z;
    }
};

struct stepper_functor
{

    odeint::controlled_runge_kutta<
        odeint::runge_kutta_cash_karp54_classic<
	    state_type, value_type, state_type, value_type,
	    odeint::vector_space_algebra, odeint::default_operations,
	    odeint::never_resizer
	    > > stepper;

//    value_type t , dt;

    stepper_functor( void ) /*: t( 0.0 ) , dt( 0.01 ) */ { }

    template <class T>
    __host__ __device__ void operator()(T s)
    {
        using namespace odeint;

        state_type    &state = thrust::get<0>(s);
        lorenz_system &sys   = thrust::get<1>(s);
        value_type &t = thrust::get<2>(s);
        value_type &dt = thrust::get<3>(s);
        

        const size_t max_attempts = 1000;

        size_t count = 0;
        while( t < t_max )
        {
            if( t_max < ( t + dt ) )
            {
                dt = t_max - t;
            }

            size_t trials = 0;
            controlled_step_result res = success;
            do
            {
                res = stepper.try_step( sys , state , t , dt );
                ++trials;
            }
            while( ( res == fail ) && ( trials < max_attempts ) );
            if( trials == max_attempts ) break;
            ++count;
            #ifndef __CUDACC__
            std::cout << t << "\t" << dt << "\t" << state << "\n";
            #endif
        }
    }
};



//---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    using namespace std;

    size_t n = argc > 1 ? atoi(argv[1]) : 1024;

    std::vector<lorenz_system> ensemble_host(n);
    value_type Rmin = 0.1 , Rmax = 50.0 , dR = ( n > 1 ) ? ( Rmax - Rmin ) / value_type( n - 1 ) : 0.0;
    for( size_t i=0 ; i<n ; ++i )
        ensemble_host[i] = lorenz_system(Rmin + dR * value_type( i ));
    thrust::device_vector<lorenz_system> ensemble = ensemble_host;

    state_type seed = {10, 10, 10};
    thrust::device_vector<state_type> x(n);
    thrust::device_vector< value_type > t(n) , dt(n);
    thrust::fill(x.begin(), x.end(), seed);
    thrust::fill(t.begin(), t.end(), 0.0);
    thrust::fill(dt.begin(), dt.end(), 0.01);


    stepper_functor step;
    thrust::for_each(
        thrust::make_zip_iterator(
            thrust::make_tuple(x.begin(), ensemble.begin(), t.begin(), dt.begin() )),
        thrust::make_zip_iterator(
            thrust::make_tuple(x.end(), ensemble.end(), t.end(), dt.end() )),
        step);

//        odeint::integrate_const(stepper, std::ref(sys[i]), X[i], double(0), t_max, dt);



    for( size_t i=0 ; i<n ; ++i )
        std::cout << ensemble_host[i].R << "\t" << x[i] << "\t" << t[i] << "\t" << dt[i] << std::endl;




    // DEBUG STUFF
    //
    // lorenz_system l( 28.0 );
    // odeint::controlled_runge_kutta<
    //     odeint::runge_kutta_cash_karp54_classic<
    //         state_type, value_type, state_type, value_type,
    //         odeint::vector_space_algebra, odeint::default_operations,
    //         odeint::never_resizer
    //         > > stepper2;
    // state_type xx = { 10.0 , 10.0 , 10.0 };
    // double t = 0.0 , dt = 0.01;
    // odeint::controlled_step_result res = stepper2.try_step( l , xx , t , dt );
    // cout << 28.0 << "\t" << t << "\t" << dt << "\t" << int( res) << "\t" << xx << endl;


    // state_type x_old = { 10.0 , 10.0 , 10.0 } , x_new;
    // state_type dxdt_old;
    // state_type x_err;
    // value_type dt = 0.01;
    // lorenz_system l( 28.0 );


    // odeint::runge_kutta_cash_karp54_classic<
    //     state_type, value_type, state_type, value_type,
    //     odeint::vector_space_algebra, odeint::default_operations,
    //     odeint::never_resizer
    //     > stepper;

    // l( x_old , dxdt_old , 0.0 );
    // stepper.do_step( l , x_old , dxdt_old , 0.0 , x_new , dt , x_err );

    // // cout << x_old << endl;
    // // cout << dxdt_old << endl;
    // // cout << x_new << endl;
    // // cout << x_err << endl;
    
    // // value_type eps_abs = 1.0e-6 , eps_rel = 1.0e-6 , a_x = 1.0 , a_dxdt = 1.0;

    // // using namespace odeint;
    // // vector_space_algebra algebra;
    // // algebra.for_each3( x_err , x_old , dxdt_old ,
    // //                    default_operations::rel_error< value_type >( eps_abs , eps_rel , a_x , a_dxdt * dt ) );

    // // value_type res = algebra.reduce( x_err , default_operations::maximum< value_type >() , 0.0 );

    // // cout << x_err << endl;
    // // cout << res << endl;

    // odeint::controlled_runge_kutta<
    //     odeint::runge_kutta_cash_karp54_classic<
    //         state_type, value_type, state_type, value_type,
    //         odeint::vector_space_algebra, odeint::default_operations,
    //         odeint::never_resizer
    //         > > stepper2;


    // double t = 0.0;
    // odeint::controlled_step_result res = stepper2.try_step( l , x_old , dxdt_old , t , x_new , dt );

    // cout << int( res ) << " " << t << " " << dt << endl;
    // cout << x_old << endl;
    // cout << dxdt_old << endl;
    // cout << x_new << endl;

}
