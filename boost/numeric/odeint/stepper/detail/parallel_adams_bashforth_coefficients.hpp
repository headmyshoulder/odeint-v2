/*
 [auto_generated]
 boost/numeric/odeint/stepper/detail/parallel_adams_bashforth_coefficients.hpp

 [begin_description]
 Definition of the coefficients for the parallel Adams-Bashforth method.
 [end_description]

 Copyright 2009-2013 Karsten Ahnert
 Copyright 2009-2013 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_PARALLEL_ADAMS_BASHFORTH_COEFFICIENTS_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_PARALLEL_ADAMS_BASHFORTH_COEFFICIENTS_HPP_INCLUDED

#include <boost/array.hpp>
#include <cmath>

using std::sqrt;

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< class Value , size_t Stages >
class parallel_adams_bashforth_coefficients ;

template< class Value >
class parallel_adams_bashforth_coefficients< Value , 2 > : public boost::array< Value , 2 >
{
public:
    parallel_adams_bashforth_coefficients( void )
        : boost::array< Value , 2 >()
    {
        (*this)[0] = static_cast< Value >(3) / static_cast< Value >(2);
        (*this)[1] = static_cast< Value >( 1 );
    }
};

template< class Value >
class parallel_adams_bashforth_coefficients< Value , 3 > : public boost::array< Value , 3 >
{
public:
    parallel_adams_bashforth_coefficients( void )
        : boost::array< Value , 3 >()
    {
        (*this)[0] = (static_cast< Value >(16) - sqrt(static_cast< Value >(6))) / static_cast< Value >(10);
        (*this)[1] = (static_cast< Value >(16) + sqrt(static_cast< Value >(6))) / static_cast< Value >(10);
        (*this)[2] = static_cast< Value >( 1 );
    }
};

template< class Value >
class parallel_adams_bashforth_coefficients< Value , 4 > : public boost::array< Value , 4 >
{
public:
    parallel_adams_bashforth_coefficients( void )
        : boost::array< Value , 4 >()
    {
        (*this)[0] = static_cast< Value >(2);
        (*this)[1] = (static_cast< Value >(15) + sqrt(static_cast< Value >(5))) / static_cast< Value >(10);
        (*this)[2] = (static_cast< Value >(15) - sqrt(static_cast< Value >(5))) / static_cast< Value >(10);
        (*this)[3] = static_cast< Value >( 1 );
    }
};

template< class Value >
class parallel_adams_bashforth_coefficients< Value , 5 > : public boost::array< Value , 5 >
{
public:
    parallel_adams_bashforth_coefficients( void )
        : boost::array< Value , 5 >()
    {
        (*this)[0] = static_cast< Value >(2);
        (*this)[1] = (static_cast< Value >(21) + sqrt(static_cast< Value >(21))) / static_cast< Value >(14);
        (*this)[2] = static_cast< Value >(3) / static_cast< Value >(2);
        (*this)[3] = (static_cast< Value >(21) - sqrt(static_cast< Value >(21))) / static_cast< Value >(14);
        (*this)[4] = static_cast< Value >( 1 );
    }
};

} // detail
} // odeint
} // numeric
} // boost



#endif // BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_PARALLEL_ADAMS_BASHFORTH_COEFFICIENTS_HPP_INCLUDED
