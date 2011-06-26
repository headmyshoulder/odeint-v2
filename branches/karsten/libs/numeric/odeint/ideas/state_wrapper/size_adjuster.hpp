#include <boost/numeric/odeint/util/is_resizeable.hpp>


template< class Stepper , class State >
bool adjust_size_by_resizeability( Stepper &stepper , const State &x , boost::true_type )
{
    return stepper.adjust_size( x );
}

template< class Stepper , class State>
bool adjust_size_by_resizeability( Stepper &stepper , const State &x , boost::false_type )
{
    return false;
}

struct always_resizer
{

    template< class Stepper , class State >
    bool adjust_size( Stepper& stepper, const State &x )
    {
        return adjust_size_by_resizeability( stepper , x , typename boost::numeric::odeint::is_resizeable< State >::type() );
        stepper.resize( x );
    }

};


struct initially_resizer
{
    bool m_initialized;

    initially_resizer(): m_initialized( false )
    { }

    template< class Stepper , class State >
    bool adjust_size( Stepper& stepper, const State &x )
    {
        if( !m_initialized )
        {
            m_initialized = true;
            return adjust_size_by_resizeability( stepper , x , typename boost::numeric::odeint::is_resizeable< State >::type() );
        }
        return false;
    }
};


struct never_resizer
{
    template< class Stepper , class State >
    void adjust_size( Stepper& stepper, const State &x )
    { }
};
