/* example for self defined vector type */

#include <vector>

#include <boost/numeric/odeint.hpp>

//[my_vector
class my_vector : public std::vector< double >
{
public:
    my_vector( const size_t N )
        : std::vector< double >( N )
    { }

    my_vector()
        : std::vector< double >()
    { }

// ...
};

// define my_vector as reizeable

namespace boost { namespace numeric { namespace odeint {

template<>
struct is_resizeable< my_vector >
{
    typedef boost::true_type type;
    static const bool value = type::value;
};

} } }
//]


typedef my_vector state_type;

void lorenz( const state_type &x , state_type &dxdt , const double t )
{
    const double sigma( 10.0 );
    const double R( 28.0 );
    const double b( 8.0 / 3.0 );

    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = -b * x[2] + x[0] * x[1];
}

using namespace boost::numeric::odeint;

int main()
{
    state_type x(3);
    x[0] = 5.0 ; x[1] = 10.0 ; x[2] = 10.0;

    // my_vector works with range_algebra as it's derived from std::vector
    // no further work is required

    integrate_const( runge_kutta4< state_type >() , lorenz , x , 0.0 , 10.0 , 0.1 );
}
