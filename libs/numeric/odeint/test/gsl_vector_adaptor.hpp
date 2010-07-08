/*
 boost header: xyz/gsl_vector_adaptor.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef GSL_VECTOR_ADAPTOR_HPP_INCLUDED
#define GSL_VECTOR_ADAPTOR_HPP_INCLUDED

#include <gsl/gsl_vector.h>

#include <boost/range.hpp>

namespace boost
{
	// ToDo define gsl_vector_iterator which increments x by stride

	template<>
    struct range_iterator< gsl_vector >
	{
		typedef double* type;
	};

	template<>
	range_iterator< gsl_vector >::type begin( gsl_vector &r )
	{
		return r.data;
	}

	template<>
	range_iterator< gsl_vector >::type end( gsl_vector &r )
	{
		return r.data + r.size;
	}



}

#endif // GSL_VECTOR_ADAPTOR_HPP_INCLUDED
