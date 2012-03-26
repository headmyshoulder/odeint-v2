#include <iostream>
#include <random>
#include <array>

#include <boost/numeric/odeint.hpp>


struct sys
{
    template< class State >
    void operator()( const State &x , State &dxdt )
    {
        dxdt[0] = -x[0];
    }
};

class stochastic_euler
{
public:

    typedef boost::numeric::odeint::stepper_tag stepper_category;



    template< class System , class State , class Time >
    void do_step( System system , State &x , Time t , Time dt ) const
    {
        State dxdt;
        system.first( x , dxdt );
        for( size_t i=0 ; i<x.size() ; ++i )
            x[i] += dt * dxdt[i] + sqrt( dt ) * system.second();
    }
};

template< class Rng , class Dist >
struct gen
{
    Rng &rng;
    Dist &dist;
    gen( Rng &rng_ , Dist &dist_ ) : rng( rng_ ) , dist( dist_ ) { }
    double operator()( void ) { return dist( rng ); }
};

template< class Rng , class Dist >
gen< Rng , Dist > make_gen( Rng &rng , Dist &dist )
{
    return gen< Rng , Dist >( rng , dist );
}


struct streaming_observer
{
    template< class State >
    void operator()( const State &x , double t ) const
    {
        std::cout << t << "\t" << x[0] << "\n";
    }
};

using namespace std;
using namespace boost::numeric::odeint;

typedef std::array< double , 1 > state_type;

int main( int argc , char **argv )
{

    mt19937 rng;
    normal_distribution<> dist( 0.0 , 1.0 );


    double dt = 0.1;
    state_type x = {{ 1.0 }};
    streaming_observer obs;
    cout << boost::is_same< streaming_observer , boost::numeric::odeint::null_observer >::value << endl;
    integrate_const( stochastic_euler() , make_pair( sys() , make_gen( rng , dist ) ) , x , 0.0 , 10.0 , dt , obs );

    return 0;
}
