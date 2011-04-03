/*
 * taylor.hpp
 *
 *  Created on: Apr 2, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_TAYLOR_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_TAYLOR_HPP_

#include <tr1/array>

// general boost includes
#include <boost/static_assert.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>

// fusion includes
#include <boost/fusion/include/is_sequence.hpp>
#include <boost/fusion/include/size.hpp>
#include <boost/fusion/include/at.hpp>

// mpl includes
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/for_each.hpp>

// proto includes
#include <boost/proto/proto.hpp>




namespace boost {
namespace numeric {
namespace odeint {


namespace taylor_adf
{

	namespace proto = boost::proto;

	template<int I> struct placeholder {};

	proto::terminal< placeholder< 0 > >::type const arg1 = {{}};
	proto::terminal< placeholder< 1 > >::type const arg2 = {{}};
	proto::terminal< placeholder< 2 > >::type const arg3 = {{}};



	template< class State , size_t Order >
	struct context_derivs : proto::callable_context< context_derivs< State , Order > const >
	{
		typedef double result_type;

		const State &m_x;
		std::tr1::array< State , Order > &m_derivs;
		size_t m_which;

		context_derivs( const State &x , std::tr1::array< State , Order > &derivs , size_t which )
		: m_x( x ) , m_derivs( derivs ) , m_which( which ) { }



		template< int I >
		double operator()( proto::tag::terminal , placeholder<I> ) const
		{
			return ( m_which == 0 ) ? m_x[I] : m_derivs[m_which-1][I];
		}

		double operator()( proto::tag::terminal , double x ) const
		{
			return ( m_which == 0 ) ? x : 0.0;
		}


		template< typename L , typename R >
		double operator ()( proto::tag::plus , L const &l , R const &r ) const
		{
			return proto::eval( l , context_derivs< State , Order >( m_x , m_derivs , m_which ) ) +
				proto::eval( r , context_derivs< State , Order >( m_x , m_derivs , m_which ) );
		}


		template< typename L , typename R >
		double operator ()( proto::tag::minus , L const &l , R const &r ) const
		{
			return proto::eval( l , context_derivs< State , Order >( m_x , m_derivs , m_which ) ) -
				proto::eval( r , context_derivs< State , Order >( m_x , m_derivs , m_which ) );
		}


		template< typename L , typename R >
		double operator ()( proto::tag::multiplies , L const &l , R const &r ) const
		{
			double tmp = 0.0;
			for( size_t k=0 ; k<=m_which ; ++k )
			{
				tmp += boost::math::binomial_coefficient< double >( m_which , k )
					* proto::eval( l , context_derivs< State , Order >( m_x , m_derivs , k ) )
					* proto::eval( r , context_derivs< State , Order >( m_x , m_derivs , m_which - k ) );
			}
			return tmp;
		}


		template< typename L , typename R >
		double operator ()( proto::tag::divides , L const &l , R const &r ) const
		{
			return 0.0;
		}

	};

	template< class System , class State , size_t Order >
	struct eval_derivs
	{
		typedef std::tr1::array< State , Order > derivs_type;

		System m_sys;
		const State &m_x;
		derivs_type &m_derivs;
		size_t m_which;
		double m_dt;

		eval_derivs( System sys , const State &x , derivs_type &derivs , size_t which , double dt )
		: m_sys( sys ) , m_x( x ) , m_derivs( derivs ) , m_which( which ) , m_dt( dt ) { }

		template< class Index >
		void operator()( Index )
		{
			m_derivs[m_which][ Index::value ] = m_dt / double( m_which + 1 ) * boost::proto::eval(
					boost::fusion::at< Index >( m_sys ) ,
					taylor_adf::context_derivs< State , Order >( m_x , m_derivs , m_which ) );
		}
	};

	template< class System , class State , size_t Order >
	eval_derivs< System , State , Order > make_eval_derivs( System sys , const State &x , std::tr1::array< State , Order > &derivs , size_t i , double dt )
	{
		return eval_derivs< System , State , Order >( sys , x , derivs , i , dt );
	}



}

template< size_t N , size_t Order >
class taylor
{
public:

	static const size_t dim = N;
	static const size_t order = Order;
	typedef double value_type;
	typedef value_type time_type;

	typedef std::tr1::array< value_type , dim > state_type;
	typedef std::tr1::array< state_type , order > derivs_type;



	template< class System >
	void do_step( System sys , state_type &x , value_type t , value_type dt )
	{

		BOOST_STATIC_ASSERT(( boost::fusion::traits::is_sequence< System >::value ));
		BOOST_STATIC_ASSERT(( size_t( boost::fusion::result_of::size< System >::value ) == dim ));

		for( size_t i=0 ; i<Order ; ++i )
		{
			boost::mpl::for_each< boost::mpl::range_c< size_t , 0 , dim > >( make_eval_derivs( sys , x , m_derivs , i , dt ) );
		}

		for( size_t i=0 ; i<N ; ++i )
		{
			for( size_t k=0 ; k<order ; ++k )
				x[i] += m_derivs[k][i];
		}
	}

	const derivs_type& get_last_derivs( void ) const
	{
		return m_derivs;
	}



private:

	derivs_type m_derivs;

};



namespace taylor_adf2
{

	namespace proto = boost::proto;

	template<int I> struct placeholder {};

	proto::terminal< placeholder< 0 > >::type const arg1 = {{}};
	proto::terminal< placeholder< 1 > >::type const arg2 = {{}};
	proto::terminal< placeholder< 2 > >::type const arg3 = {{}};



	template< class State , size_t Order >
	struct context_derivs : proto::callable_context< context_derivs< State , Order > const >
	{
		typedef double result_type;

		const State &m_x;
		std::tr1::array< State , Order > &m_derivs;
		size_t m_which;

		context_derivs( const State &x , std::tr1::array< State , Order > &derivs , size_t which )
		: m_x( x ) , m_derivs( derivs ) , m_which( which ) { }



		template< int I >
		double operator()( proto::tag::terminal , placeholder<I> ) const
		{
			return ( m_which == 0 ) ? m_x[I] : m_derivs[m_which-1][I];
		}

		double operator()( proto::tag::terminal , double x ) const
		{
			return ( m_which == 0 ) ? x : 0.0;
		}


		template< typename L , typename R >
		double operator ()( proto::tag::plus , L const &l , R const &r ) const
		{
			return proto::eval( l , context_derivs< State , Order >( m_x , m_derivs , m_which ) ) +
				proto::eval( r , context_derivs< State , Order >( m_x , m_derivs , m_which ) );
		}


		template< typename L , typename R >
		double operator ()( proto::tag::minus , L const &l , R const &r ) const
		{
			return proto::eval( l , context_derivs< State , Order >( m_x , m_derivs , m_which ) ) -
				proto::eval( r , context_derivs< State , Order >( m_x , m_derivs , m_which ) );
		}


		template< typename L , typename R >
		double operator ()( proto::tag::multiplies , L const &l , R const &r ) const
		{
			double tmp = 0.0;
			for( size_t k=0 ; k<=m_which ; ++k )
			{
				tmp += boost::math::binomial_coefficient< double >( m_which , k )
					* proto::eval( l , context_derivs< State , Order >( m_x , m_derivs , k ) )
					* proto::eval( r , context_derivs< State , Order >( m_x , m_derivs , m_which - k ) );
			}
			return tmp;
		}


		template< typename L , typename R >
		double operator ()( proto::tag::divides , L const &l , R const &r ) const
		{
			return 0.0;
		}

	};

	template< class System , class State , size_t Order >
	struct eval_derivs
	{
		typedef std::tr1::array< State , Order > derivs_type;

		System m_sys;
		const State &m_x;
		derivs_type &m_derivs;
		size_t m_which;

		eval_derivs( System sys , const State &x , derivs_type &derivs , size_t which )
		: m_sys( sys ) , m_x( x ) , m_derivs( derivs ) , m_which( which ) { }

		template< class Index >
		void operator()( Index )
		{
			m_derivs[m_which][ Index::value ] = boost::proto::eval(
					boost::fusion::at< Index >( m_sys ) ,
					taylor_adf::context_derivs< State , Order >( m_x , m_derivs , m_which ) );
		}
	};

	template< class System , class State , size_t Order >
	eval_derivs< System , State , Order > make_eval_derivs( System sys , const State &x , std::tr1::array< State , Order > &derivs , size_t i )
	{
		return eval_derivs< System , State , Order >( sys , x , derivs , i );
	}



}




} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_STEPPER_TAYLOR_HPP_ */
