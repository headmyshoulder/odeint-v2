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

template< class Value >
class parallel_adams_bashforth_coefficients< Value , 6 > : public boost::array< Value , 6 >
{
public:
    parallel_adams_bashforth_coefficients( void )
        : boost::array< Value , 6 >()
    {
        (*this)[0] = static_cast< Value >(2);
        // x = 1 + 1/42 (21+sqrt(21 (7+2 sqrt(7))))
        (*this)[1] = static_cast< Value >(1)/static_cast< Value >(42) * ( static_cast< Value >(63)+sqrt(static_cast< Value >(21)*(static_cast< Value >(7)+static_cast< Value >(2)*sqrt(static_cast< Value >(7)))));
        // x = 1/42 (21+sqrt(21 (7-2 sqrt(7))))
        (*this)[2] = static_cast< Value >(1)/static_cast< Value >(42) * ( static_cast< Value >(63)+sqrt(static_cast< Value >(21)*(static_cast< Value >(7)-static_cast< Value >(2)*sqrt(static_cast< Value >(7)))));
        (*this)[3] = static_cast< Value >(3) - (*this)[1];
        (*this)[4] = static_cast< Value >(3) - (*this)[2];
        (*this)[5] = static_cast< Value >(1);
    }
};

template< class Value >
class parallel_adams_bashforth_coefficients< Value , 7 > : public boost::array< Value , 7 >
{
public:
    parallel_adams_bashforth_coefficients( void )
        : boost::array< Value , 7 >()
    {
        (*this)[0] = static_cast<Value>(2);
        // x = 1 + 1/66 (33+sqrt(495+66 sqrt(15)))
        (*this)[1] = static_cast<Value>(1)/static_cast<Value>(66) * (static_cast<Value>(99) + sqrt(static_cast<Value>(495)+static_cast<Value>(66)*sqrt(static_cast<Value>(15))));
        // x = 1 + 1/66 (33+sqrt(495-66 sqrt(15)))
        (*this)[2] = static_cast<Value>(1)/static_cast<Value>(66) * (static_cast<Value>(99) + sqrt(static_cast<Value>(495)-static_cast<Value>(66)*sqrt(static_cast<Value>(15))));
        (*this)[3] = static_cast< Value >(3)/static_cast<Value>(2);
        (*this)[4] = static_cast< Value >(3) - (*this)[1];
        (*this)[5] = static_cast< Value >(3) - (*this)[2];
        (*this)[6] = static_cast< Value >(1);
    }
};

} // detail
} // odeint
} // numeric
} // boost



#endif // BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_PARALLEL_ADAMS_BASHFORTH_COEFFICIENTS_HPP_INCLUDED
