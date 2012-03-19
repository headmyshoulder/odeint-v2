/*
 * units_helper.hpp
 *
 *  Created on: Mar 19, 2012
 *      Author: mario
 */

#ifndef UNITS_HELPER_HPP_
#define UNITS_HELPER_HPP_

#include <boost/units/quantity.hpp>

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< class T >
struct unit_value_type {
    typedef T type;
};


template< class Unit , class Y >
struct unit_value_type< boost::units::quantity<Unit , Y > > {
    typedef Y type;
};


}
}
}
}

#endif /* UNITS_HELPER_HPP_ */
