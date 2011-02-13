/*
 * is_pair.hpp
 *
 *  Created on: Feb 12, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_IS_PAIR_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_IS_PAIR_HPP_

#include <boost/mpl/bool.hpp>
#include <utility>


namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< class T >
struct is_pair : public boost::mpl::false_
{
};

template< class T1 , class T2 >
struct is_pair< std::pair< T1 , T2 > > : public boost::mpl::true_
{
};

} // namespace detail
} // namespace odeint
} // namespace numeric
} // namespace boost

#endif /* BOOST_NUMERIC_ODEINT_STEPPER_DETAIL_IS_PAIR_HPP_ */
