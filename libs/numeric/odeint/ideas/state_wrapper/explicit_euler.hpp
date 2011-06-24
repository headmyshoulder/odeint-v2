#include <iostream>

#include <boost/range.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/util/size_adjuster.hpp>

template< class V >
struct state_wrapper
{
    typedef typename V::value_type value_type;

    V m_v;

    state_wrapper() : m_v() { }

};


template< typename StateType >
class explicit_euler {

public:

    typedef double value_type;
    typedef double time_type;
    typedef StateType state_type;
    typedef state_wrapper< state_type > wrapped_state_type;

    explicit_euler() : m_dxdt()
    {
        m_size_adjuster.register_state( 0 , m_dxdt.m_v );
    };

    template< class System , class StateInOut >
    void do_step( System system , StateInOut &inout , const time_type &t , const time_type &dt )
    {
        m_size_adjuster.adjust_size( inout );

        system( inout , m_dxdt.m_v , t );

        boost::numeric::odeint::range_algebra::for_each3( inout , inout , m_dxdt.m_v , typename boost::numeric::odeint::default_operations::template scale_sum2< value_type , time_type >( 1.0 , dt ) );
    }

private:
    boost::numeric::odeint::size_adjuster< state_type , 1 > m_size_adjuster;
    wrapped_state_type m_dxdt;
};
