/*
 * generic_rk_operations.hpp
 *
 *  Created on: May 30, 2011
 *      Author: mario
 */

#ifndef GENERIC_RK_OPERATIONS_HPP_
#define GENERIC_RK_OPERATIONS_HPP_

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< size_t StageNumber , class Operations , class Fac >
struct generic_rk_scale_sum;

template< class Operations , class Fac >
struct generic_rk_scale_sum< 1 , Operations , Fac > : public Operations::template scale_sum2< Fac >
{
    generic_rk_scale_sum( const boost::array<Fac,1> &a , const Fac &dt ) : Operations::template scale_sum2< Fac >( 1.0 , a[0]*dt )
        { }

    typedef void result_type;
};


template< class Operations , class Fac >
struct generic_rk_scale_sum< 2 , Operations , Fac > : public Operations::template scale_sum3< Fac >
{
    generic_rk_scale_sum( const boost::array<Fac,2> &a , const Fac &dt )
            : Operations::template scale_sum3< Fac >( 1.0 , a[0]*dt , a[1]*dt )
        { }

    typedef void result_type;
};

template< class Operations , class Fac >
struct generic_rk_scale_sum< 3 , Operations , Fac > : public Operations::template scale_sum4< Fac >
{
    generic_rk_scale_sum( const boost::array<Fac,3> &a , const Fac &dt )
            : Operations::template scale_sum4< Fac >( 1.0 , a[0]*dt , a[1]*dt , a[2]*dt )
        { }

    typedef void result_type;
};

template< class Operations , class Fac >
struct generic_rk_scale_sum< 4 , Operations , Fac > : public Operations::template scale_sum5< Fac >
{
    generic_rk_scale_sum( const boost::array<Fac,4> &a , const Fac &dt )
            : Operations::template scale_sum5< Fac >( 1.0 , a[0]*dt , a[1]*dt , a[2]*dt , a[3]*dt )
        { }

    typedef void result_type;
};

template< class Operations , class Fac >
struct generic_rk_scale_sum< 5 , Operations , Fac > : public Operations::template scale_sum6< Fac >
{
    generic_rk_scale_sum( const boost::array<Fac,5> &a , const Fac &dt )
            : Operations::template scale_sum6< Fac >( 1.0 , a[0]*dt , a[1]*dt , a[2]*dt , a[3]*dt , a[4]*dt )
        { }

    typedef void result_type;
};

template< class Operations , class Fac >
struct generic_rk_scale_sum< 6 , Operations , Fac > : public Operations::template scale_sum7< Fac >
{
    generic_rk_scale_sum( const boost::array<Fac,6> &a , const Fac &dt )
            : Operations::template scale_sum7< Fac >( 1.0 , a[0]*dt , a[1]*dt , a[2]*dt , a[3]*dt , a[4]*dt , a[5]*dt )
        { }

    typedef void result_type;
};

// for error estimates
template< size_t StageNumber , class Operations , class Fac >
struct generic_rk_scale_sum_err;

template< class Operations , class Fac >
struct generic_rk_scale_sum_err< 1 , Operations , Fac > : public Operations::template scale_sum1< Fac >
{
    generic_rk_scale_sum_err( const boost::array<Fac,1> &a , const Fac &dt ) : Operations::template scale_sum1< Fac >( a[0]*dt )
        { }

    typedef void result_type;
};


template< class Operations , class Fac >
struct generic_rk_scale_sum_err< 2 , Operations , Fac > : public Operations::template scale_sum2< Fac >
{
    generic_rk_scale_sum_err( const boost::array<Fac,2> &a , const Fac &dt )
            : Operations::template scale_sum2< Fac >( a[0]*dt , a[1]*dt )
        { }

    typedef void result_type;
};

template< class Operations , class Fac >
struct generic_rk_scale_sum_err< 3 , Operations , Fac > : public Operations::template scale_sum3< Fac >
{
    generic_rk_scale_sum_err( const boost::array<Fac,3> &a , const Fac &dt )
            : Operations::template scale_sum3< Fac >( a[0]*dt , a[1]*dt , a[2]*dt )
        { }

    typedef void result_type;
};

template< class Operations , class Fac >
struct generic_rk_scale_sum_err< 4 , Operations , Fac > : public Operations::template scale_sum4< Fac >
{
    generic_rk_scale_sum_err( const boost::array<Fac,4> &a , const Fac &dt )
            : Operations::template scale_sum4< Fac >( a[0]*dt , a[1]*dt , a[2]*dt , a[3]*dt )
        { }

    typedef void result_type;
};

template< class Operations , class Fac >
struct generic_rk_scale_sum_err< 5 , Operations , Fac > : public Operations::template scale_sum5< Fac >
{
    generic_rk_scale_sum_err( const boost::array<Fac,5> &a , const Fac &dt )
            : Operations::template scale_sum5< Fac >( a[0]*dt , a[1]*dt , a[2]*dt , a[3]*dt , a[4]*dt )
        { }

    typedef void result_type;
};


template< class Operations , class Fac >
struct generic_rk_scale_sum_err< 6 , Operations , Fac > : public Operations::template scale_sum6< Fac >
{
    generic_rk_scale_sum_err( const boost::array<Fac,6> &a , const Fac &dt )
            : Operations::template scale_sum6< Fac >( a[0]*dt , a[1]*dt , a[2]*dt , a[3]*dt , a[4]*dt , a[5]*dt )
        { }

    typedef void result_type;
};

}
}
}
}


#endif /* GENERIC_RK_OPERATIONS_HPP_ */
