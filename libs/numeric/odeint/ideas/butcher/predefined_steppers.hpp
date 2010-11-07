/*
 * predefined_steppers.hpp
 *
 *  Created on: Nov 7, 2010
 *      Author: karsten
 */

#ifndef PREDEFINED_STEPPERS_HPP_
#define PREDEFINED_STEPPERS_HPP_

#include <boost/mpl/int.hpp>
#include <boost/ratio.hpp>

#include "explicit_runge_kutta.hpp"

typedef mpl::int_< 0 > null;
typedef mpl::int_< 1 > one;
typedef boost::ratio< 1 , 2 > one_half;
typedef boost::ratio< 1 , 3 > one_third;
typedef boost::ratio< 1 , 6 > one_sixth;



/*
 * euler :
 * 0 |
 *   | 1
 */

typedef mpl::vector< null > euler_c;
typedef mpl::vector<> euler_a;
typedef mpl::vector< one > euler_b;

template< class state_type >
struct mpl_euler_stepper : public explicit_runge_kutta< state_type , euler_a , euler_c , euler_b > { };


/*
 * midpoint :
 * 0   |
 * 1/2 | 1/2
 * -----------
 *     | 0   1
 */

typedef mpl::vector< null , one_half > midpoint_c;
typedef mpl::vector< mpl::vector< one_half > > midpoint_a;
typedef mpl::vector< null , one > midpoint_b;

template< class state_type >
struct mpl_midpoint_stepper : public explicit_runge_kutta< state_type , midpoint_a , midpoint_c , midpoint_b > { };




/*
 * classical Runge Kutta
 * 0   |
 * 1/2 | 1/2
 * 1/2 | 0   1/2
 * 1   | 0   0   1
 * ------------------
 *     | 1/6 1/3 1/6 1/6
 *
 */
typedef mpl::vector< null , one_half , one_half , one > rk4_c;
typedef mpl::vector<
			mpl::vector< one_half > ,
			mpl::vector< null , one_half > ,
			mpl::vector< null , null , one >
		> rk4_a;
typedef mpl::vector< one_sixth , one_third , one_third , one_sixth > rk4_b;

template< class state_type >
struct mpl_rk4_stepper : public explicit_runge_kutta< state_type , rk4_a , rk4_c , rk4_b > { };

#endif /* PREDEFINED_STEPPERS_HPP_ */
