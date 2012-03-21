/*
 [auto_generated]
 boost/numeric/odeint/util/number_of_elements.hpp

 [begin_description]
 Generic number of elements (size) function.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_UTIL_NUMBER_OF_ELEMENTS_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_UTIL_NUMBER_OF_ELEMENTS_HPP_INCLUDED



namespace boost {
namespace numeric {
namespace odeint {

template< class T >
struct number_of_elements_impl
{
    typedef typename T::size_type size_type;
    typedef size_type result_type;

    static result_type number_of_elements( const T &x )
    {
        return x.size();
    }


};


template< class T >
typename number_of_elements_impl< T >::result_type number_of_elements( const T &x )
{
    return number_of_elements_impl< T >::number_of_elements( x );
}


} // namespace odeint
} // namespace numeric
} // namespace boost



#endif // BOOST_NUMERIC_ODEINT_UTIL_NUMBER_OF_ELEMENTS_HPP_INCLUDED
