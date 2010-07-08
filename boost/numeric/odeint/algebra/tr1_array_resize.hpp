/*
 boost header: BOOST_NUMERIC_ODEINT/tr1_array_resize.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_TR1_ARRAY_RESIZE_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_TR1_ARRAY_RESIZE_HPP_INCLUDED


namespace boost {
namespace numeric {
namespace odeint {

template< class T , size_t N >
struct is_resizeable< std::tr1::array< T , N > >
{
	struct type : public boost::false_type { };
	const static bool value = type::value;
};


} // odeint
} // numeric
} // boost


#endif //BOOST_BOOST_NUMERIC_ODEINT_TR1_ARRAY_RESIZE_HPP_INCLUDED
