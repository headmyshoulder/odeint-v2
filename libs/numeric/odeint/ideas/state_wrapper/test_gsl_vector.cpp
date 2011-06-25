#include <iostream>
#include <gsl/gsl_vector.h>

#include "explicit_euler.hpp"
#include <boost/numeric/odeint/external/gsl/gsl_vector_adaptor.hpp>

using namespace std;

typedef gsl_vector state_type;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

void lorenz( const state_type x , state_type dxdt , double t )
{
    gsl_vector_set( dxdt , 0 , sigma * ( gsl_vector_get(x , 1 ) - gsl_vector_get( x , 0 ) ) );
    gsl_vector_set( dxdt , 1 , R * gsl_vector_get( x , 0 ) - gsl_vector_get( x , 1 ) - gsl_vector_get( x , 0 ) * gsl_vector_get( x , 2) );
    gsl_vector_set( dxdt , 2 , gsl_vector_get( x , 0 ) * gsl_vector_get( x , 1 ) - b * gsl_vector_get( x , 2) );
}

template<>
struct state_wrapper< state_type >
{
    typedef double value_type;

    state_type m_v;

    state_wrapper( )
    {
        m_v->owner = 0;
        m_v->size = 0;
        m_v->stride = 0;
        m_v->data = 0;
        m_v->block = 0;
    }

};

int main() {

    explicit_euler< state_type > euler;

    state_type x = gsl_vector_alloc( 3 );
    gsl_vector_set( x , 0 , 1.0);
    gsl_vector_set( x , 1 , 1.0);
    gsl_vector_set( x , 2 , 2.0);

    euler.do_step( lorenz , x , 0.0 , 0.1 );

    cout << gsl_vector_get( x , 0 ) << "  " << gsl_vector_get( x , 1 ) << "  " << gsl_vector_get( x , 2 ) << endl;

    gsl_vector_free( x );
}
