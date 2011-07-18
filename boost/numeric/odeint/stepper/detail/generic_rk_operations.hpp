/*
 [auto_generated]
 boost/numeric/odeint/stepper/detail/generic_rk_operations.hpp

 [begin_description]
 Operations caller for the generic Runge Kutta method.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_GENERIC_RK_OPERATIONS_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_GENERIC_RK_OPERATIONS_HPP_INCLUDED


namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< size_t StageNumber , class Operations , class Fac , class Time >
struct generic_rk_scale_sum;

template< class Operations , class Fac , class Time >
struct generic_rk_scale_sum< 1 , Operations , Fac , Time > : public Operations::template scale_sum2< Fac , Time >
{
    generic_rk_scale_sum( const boost::array<Fac,1> &a , const Time &dt ) : Operations::template scale_sum2< Fac , Time >( 1.0 , a[0]*dt )
                { }

    typedef void result_type;
};


template< class Operations , class Fac , class Time >
struct generic_rk_scale_sum< 2 , Operations , Fac , Time > : public Operations::template scale_sum3< Fac , Time >
{
    generic_rk_scale_sum( const boost::array<Fac,2> &a , const Time &dt )
                    : Operations::template scale_sum3< Fac , Time >( 1.0 , a[0]*dt , a[1]*dt )
                      { }

    typedef void result_type;
};

template< class Operations , class Fac , class Time >
struct generic_rk_scale_sum< 3 , Operations , Fac , Time > : public Operations::template scale_sum4< Fac , Time >
{
    generic_rk_scale_sum( const boost::array<Fac,3> &a , const Time &dt )
                    : Operations::template scale_sum4< Fac , Time >( 1.0 , a[0]*dt , a[1]*dt , a[2]*dt )
                      { }

    typedef void result_type;
};

template< class Operations , class Fac , class Time >
struct generic_rk_scale_sum< 4 , Operations , Fac , Time > : public Operations::template scale_sum5< Fac , Time >
{
    generic_rk_scale_sum( const boost::array<Fac,4> &a , const Time &dt )
                    : Operations::template scale_sum5< Fac , Time >( 1.0 , a[0]*dt , a[1]*dt , a[2]*dt , a[3]*dt )
                      { }

    typedef void result_type;
};

template< class Operations , class Fac , class Time >
struct generic_rk_scale_sum< 5 , Operations , Fac , Time > : public Operations::template scale_sum6< Fac , Time >
{
    generic_rk_scale_sum( const boost::array<Fac,5> &a , const Time &dt )
                    : Operations::template scale_sum6< Fac , Time >( 1.0 , a[0]*dt , a[1]*dt , a[2]*dt , a[3]*dt , a[4]*dt )
                      { }

    typedef void result_type;
};

template< class Operations , class Fac , class Time >
struct generic_rk_scale_sum< 6 , Operations , Fac , Time > : public Operations::template scale_sum7< Fac , Time >
{
    generic_rk_scale_sum( const boost::array<Fac,6> &a , const Time &dt )
                    : Operations::template scale_sum7< Fac , Time >( 1.0 , a[0]*dt , a[1]*dt , a[2]*dt , a[3]*dt , a[4]*dt , a[5]*dt )
                      { }

    typedef void result_type;
};

// for error estimates
template< size_t StageNumber , class Operations , class Fac , class Time >
struct generic_rk_scale_sum_err;

template< class Operations , class Fac , class Time >
struct generic_rk_scale_sum_err< 1 , Operations , Fac , Time > : public Operations::template scale_sum1< Fac , Time >
{
    generic_rk_scale_sum_err( const boost::array<Fac,1> &a , const Time &dt ) : Operations::template scale_sum1< Fac , Time >( a[0]*dt )
                { }

    typedef void result_type;
};


template< class Operations , class Fac , class Time >
struct generic_rk_scale_sum_err< 2 , Operations , Fac , Time > : public Operations::template scale_sum2< Fac , Time >
{
    generic_rk_scale_sum_err( const boost::array<Fac,2> &a , const Time &dt )
                    : Operations::template scale_sum2< Fac , Time >( a[0]*dt , a[1]*dt )
                      { }

    typedef void result_type;
};

template< class Operations , class Fac , class Time >
struct generic_rk_scale_sum_err< 3 , Operations , Fac , Time > : public Operations::template scale_sum3< Fac , Time >
{
    generic_rk_scale_sum_err( const boost::array<Fac,3> &a , const Time &dt )
                    : Operations::template scale_sum3< Fac , Time >( a[0]*dt , a[1]*dt , a[2]*dt )
                      { }

    typedef void result_type;
};

template< class Operations , class Fac , class Time >
struct generic_rk_scale_sum_err< 4 , Operations , Fac , Time > : public Operations::template scale_sum4< Fac , Time >
{
    generic_rk_scale_sum_err( const boost::array<Fac,4> &a , const Time &dt )
                    : Operations::template scale_sum4< Fac , Time >( a[0]*dt , a[1]*dt , a[2]*dt , a[3]*dt )
                      { }

    typedef void result_type;
};

template< class Operations , class Fac , class Time >
struct generic_rk_scale_sum_err< 5 , Operations , Fac , Time > : public Operations::template scale_sum5< Fac , Time >
{
    generic_rk_scale_sum_err( const boost::array<Fac,5> &a , const Time &dt )
                    : Operations::template scale_sum5< Fac , Time >( a[0]*dt , a[1]*dt , a[2]*dt , a[3]*dt , a[4]*dt )
                      { }

    typedef void result_type;
};


template< class Operations , class Fac , class Time >
struct generic_rk_scale_sum_err< 6 , Operations , Fac , Time > : public Operations::template scale_sum6< Fac , Time >
{
    generic_rk_scale_sum_err( const boost::array<Fac,6> &a , const Time &dt )
                    : Operations::template scale_sum6< Fac , Time >( a[0]*dt , a[1]*dt , a[2]*dt , a[3]*dt , a[4]*dt , a[5]*dt )
                      { }

    typedef void result_type;
};

}
}
}
}


#endif // BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_GENERIC_RK_OPERATIONS_HPP_INCLUDED
