/*
 * explicit_runge_kutta.hpp
 *
 *  Created on: Nov 5, 2010
 *      Author: karsten
 */

#ifndef EXPLICIT_RUNGE_KUTTA_HPP_
#define EXPLICIT_RUNGE_KUTTA_HPP_

#include <boost/static_assert.hpp>
#include <boost/ratio.hpp>
#include <boost/type_traits/is_same.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/zip_view.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/at.hpp>

#include "convert_value.hpp"
#include "diagnostic.hpp"
#include "algebra.hpp"

namespace mpl = boost::mpl;


/*
 * c_1 |
 * c_2 | a_{2,1}
 * c_3 | a_{3,1} a_{3,2}
 * ...
 * c_s | a_{s,1} a_{s,2} ... a_{s,s-1}
 * -----------------------------------
 *     | b_1     b_2         b_{s-1}    b_s
 *
 * ButcherA : the coefficients a_ij
 *
 * ButcherC : the time coefficients c_i
 *
 * ButcherC : the sum coefficients b_j
 */
template< class State , class ButcherA , class ButcherC , class ButcherB >
class explicit_runge_kutta
{
public:

	typedef ButcherA butcher_a;
	typedef ButcherB butcher_b;
	typedef ButcherC butcher_c;

	static const size_t dim = mpl::size< ButcherC >::value;

	typedef mpl::range_c< size_t , 0 , dim > indices;
	typedef typename mpl::push_back< butcher_a , butcher_b >::type butcher_right;
	typedef mpl::zip_view< mpl::vector< indices , butcher_c , butcher_right > > butcher_tableau;


	typedef double time_type;
	typedef State state_type;


	/*
	 * t , dt , system, x, k[]
	 */
	template< class System >
	struct calculate_stage
	{
		System &system;
		state_type &x , &x_tmp;
		state_type *k_vector;
		const double t;
		const double dt;

		calculate_stage( System &_system , state_type &_x , state_type &_x_tmp , state_type *_k_vector , const double _t , const double _dt )
		: system( _system ) , x( _x ) , x_tmp( _x_tmp ) , k_vector( _k_vector ) , t( _t ) , dt( _dt )
		{}

		/*
		 * for( size_t j=0 ; j<dim ; ++j )
		 * {
		 * 	system( x[j] , k[j] , c[j] )
		 *  x[j+1] = sum( i , a[j][i] , k[j] )
		 * }
		 */

		template< class Stage > void operator()( Stage v )
		{
			typedef typename mpl::at< Stage , mpl::int_< 0 > >::type index_type;
			typedef typename mpl::at< Stage , mpl::int_< 1 > >::type time_value_type;
			typedef typename mpl::at< Stage , mpl::int_< 2 > >::type coef_type;

			const static size_t index = index_type::value;
			const static double time_value = convert_value< time_value_type >::get_value();

			BOOST_STATIC_ASSERT(( ( mpl::size< coef_type >::value - 1 ) == index ));

//			cout << index << " " << time_value << " " << endl;

			if( index == 0 ) system( x , k_vector[index] , t + time_value * dt );
			else system( x_tmp , k_vector[index] , t + time_value * dt );

			if( index == ( dim - 1 ) )
				algebra< state_type , coef_type , mpl::int_< index > >::do_step( x , x , k_vector , dt );
			else
				algebra< state_type , coef_type , mpl::int_< index > >::do_step( x_tmp , x , k_vector , dt );
		};
	};


	explicit_runge_kutta( void )
	{
	}

	template< class System >
	void do_step( System &system , state_type &x , double t , const double dt )
	{
		mpl::for_each< butcher_tableau >( calculate_stage< System >( system , x , m_x_tmp , m_k_vector , t , dt ) );
	}



	void print_butcher_a( void )
	{
		mpl::for_each< butcher_a >( print_vector() );
	}

	void print_butcher_b( void )
	{
		mpl::for_each< butcher_b >( print_value() );
	}

	void print_butcher_c( void )
	{
		mpl::for_each< butcher_c >( print_value() );
	}

	void print_tableau( void )
	{
		mpl::for_each< butcher_tableau >( print_tableau_vector() );
	}



private:

	state_type m_k_vector[dim];
	state_type m_x_tmp;

	BOOST_STATIC_ASSERT(( dim > 0 ));
	BOOST_STATIC_ASSERT(( dim == size_t( mpl::size< butcher_b >::value ) ));
	BOOST_STATIC_ASSERT(( ( dim - 1 ) == mpl::size< butcher_a >::value ));
};


#endif /* EXPLICIT_RUNGE_KUTTA_HPP_ */
