#include <iostream>

#include <boost/range.hpp>

#include <boost/numeric/odeint/algebra/detail/for_each.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>

template< class V >
struct state_wrapper
{
    typedef typename V::value_type value_type;

    V &m_v;
    V m_v_;

    state_wrapper( V &v ) : m_v( v ) { }

    state_wrapper() : m_v( m_v_ ) { }

    typename boost::range_iterator<V>::type begin()
    { return boost::begin( m_v ); }

    typename boost::range_iterator<V>::type end()
    { return boost::end( m_v ); }

    void resize( V &v )
    {
        m_v.resize( boost::size( v ) );
    }

    bool same_size( V &v )
    { return ( boost::size( m_v ) == boost::size( v ) ); }
};


struct range_algebra
{
    template< class S1 , class S2 , class S3 , class Op >
    static void for_each3( S1 &s1 , S2 &s2 , S3 &s3 , Op op )
    {
        boost::numeric::odeint::detail::for_each3( s1.begin() , s1.end() , s2.begin() , s3.begin() , op );
    }
};

template< typename StateType >
class explicit_euler {

public:

    typedef double value_type;
    typedef double time_type;
    typedef StateType state_type;
    typedef state_wrapper< state_type > wrapped_state_type;

    explicit_euler() { };

    template< class System , class StateInOut >
    void do_step( System system , StateInOut &inout , const time_type &t , const time_type &dt )
    {
        state_wrapper< StateInOut > x( inout );

        if( !m_dxdt.same_size( inout ) )
            m_dxdt.resize( inout );

        system( inout , m_dxdt.m_v , t );

        range_algebra::for_each3( x , x , m_dxdt , typename boost::numeric::odeint::default_operations::template scale_sum2< value_type , time_type >( 1.0 , dt ) );
    }

private:
    wrapped_state_type m_dxdt;
};
