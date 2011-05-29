#include <iostream>
#include <gsl/gsl_vector.h>

#include "explicit_euler.hpp"

using namespace std;

typedef gsl_vector *state_type;

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

    state_wrapper( state_type &v ) : m_v( v ) {
        cout << m_v->size << endl;
    }

    state_wrapper( )
    {
        m_v->owner = 0;
        m_v->size = 0;
        m_v->stride = 0;
        m_v->data = 0;
        m_v->block = 0;
    }

    double* begin()
    { return m_v->data; }

    double* end()
    { return m_v->data + m_v->size; }

    void resize( state_type &v )
    {
        cout << v->size << " " << m_v->owner << " " << v->owner << endl;

        if( m_v->owner != 0 )
        {
            gsl_block_free( m_v->block );
        }
        m_v->size = 0;

        cout << v->size << endl;

        if( v->size == 0 ) return;

        gsl_block *block = gsl_block_alloc( v->size );
        if( block == 0 ) throw std::bad_alloc( );

        m_v->data = block->data ;
        m_v->size = v->size;
        m_v->stride = 1;
        m_v->block = block;
        m_v->owner = 1;
    }

    bool same_size( state_type &v )
    { return ( m_v->size == v->size ); }
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
