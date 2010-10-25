/*
 boost header: boost/numeric/odeint/dense_output_dopri5.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_BOOST_NUMERIC_ODEINT_DENSE_OUTPUT_DOPRI5_HPP_INCLUDED
#define BOOST_BOOST_NUMERIC_ODEINT_DENSE_OUTPUT_DOPRI5_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/dense_output_dopri5.hpp>

namespace boost {
namespace numeric {
namespace odeint {


template< class ControlledStepper >
class dense_output_dopri5
{
private:

	typedef ControlledStepper stepper_type;

};

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif //BOOST_BOOST_NUMERIC_ODEINT_DENSE_OUTPUT_DOPRI5_HPP_INCLUDED
