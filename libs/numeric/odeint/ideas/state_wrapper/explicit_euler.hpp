#include <iostream>

#include <boost/range.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/util/resize.hpp>

template< class V >
struct state_wrapper
{
    typedef typename V::value_type value_type;

    V m_v;

    state_wrapper() : m_v() { }

};


template< typename StateType , class Resizer >
class explicit_euler {

public:

    typedef double value_type;
    typedef double time_type;
    typedef StateType state_type;
    typedef state_wrapper< state_type > wrapped_state_type;
    typedef Resizer resizer_type;

    explicit_euler() : m_dxdt() , m_resizer()
    { }

    template< class System , class StateInOut >
    void do_step( System system , StateInOut &inout , const time_type &t , const time_type &dt )
    {
        m_resizer.adjust_size( *this , inout );

        system( inout , m_dxdt.m_v , t );

        boost::numeric::odeint::range_algebra::for_each3( inout , inout , m_dxdt.m_v , typename boost::numeric::odeint::default_operations::template scale_sum2< value_type , time_type >( 1.0 , dt ) );
    }

    template< class State >
    bool adjust_size( const State &x )
    {
        if( boost::numeric::odeint::same_size( x , m_dxdt.m_v ) )
        {
            boost::numeric::odeint::resize( x , m_dxdt.m_v );
            return true;
        } else
            return false;
    }

private:
    wrapped_state_type m_dxdt;
    resizer_type m_resizer;
};
