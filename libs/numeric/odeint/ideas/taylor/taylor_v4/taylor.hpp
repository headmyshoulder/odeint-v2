/*
 * taylor.hpp
 *
 *  Created on: Apr 2, 2011
 *      Author: karsten
 */

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_TAYLOR_HPP_
#define BOOST_NUMERIC_ODEINT_STEPPER_TAYLOR_HPP_

#include <ostream>
#include <cmath>
#include <tr1/array>

#include <iostream>
using namespace std;
#define tab "\t"

// general boost includes
#include <boost/static_assert.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>

// fusion includes
#include <boost/fusion/include/is_sequence.hpp>
#include <boost/fusion/include/size.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/fusion/include/transform.hpp>

// mpl includes
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/for_each.hpp>

// proto includes
#include <boost/proto/proto.hpp>

/*
 * OK # dt entfernen aus context
 * OK # do_step zu try_step
 * OK # Skalierung
 *   OK -- dt_fac # Namen für Skalierungsfaktor finden
 *   OK # Faktor einschleusen in context
 *   OK # Nach jedem Durchgang einer Ordnung Skalierungsfaktor neu berechnen und schon
 *     berechnete Ableitungen skalieren
 * # Fehler als Member parameter
 * # Fehler anhand der hoechsten Ableitung berechnen
 * # Schrittweite bestimmen
 * # Schritt ausfuehren
 */


namespace boost {
namespace numeric {
namespace odeint {


namespace taylor_adf
{

	namespace proto = boost::proto;
	namespace fusion = boost::fusion;
	namespace proto = boost::proto;

	template< typename I > struct placeholder : I {};

	proto::terminal< placeholder< mpl::size_t< 0 > > >::type const arg1 = {};
	proto::terminal< placeholder< mpl::size_t< 1 > > >::type const arg2 = {};
	proto::terminal< placeholder< mpl::size_t< 2 > > >::type const arg3 = {};

	template< typename I >
	std::ostream& operator<<( std::ostream &s , const placeholder< I > &p )
	{
		s << "placeholder<" << I::value << ">";
		return s;
	}


	/*
	 * eventuel array< double , n > dt mit potenzen anlegen
	 */
	template< class State , size_t MaxOrder >
	struct taylor_context
	{
		typedef State state_type;
		typedef std::tr1::array< state_type , MaxOrder > deriv_type;

		size_t which;
		const state_type &x;
		deriv_type &derivs;
		double &dt_fac;

		taylor_context( size_t _which , const state_type &_x , deriv_type &_derivs , double &_dt_fac )
		: which( _which ) , x( _x ) , derivs( _derivs ) , dt_fac( _dt_fac ) { }
	};


	template< typename Grammar >
	struct plus_transform : proto::transform< plus_transform< Grammar > >
	{
		template< typename Expr , typename State , typename Data >
		struct impl : proto::transform_impl< Expr , State , Data >
		{
			typedef double result_type;

			result_type operator ()(
					typename impl::expr_param expr ,
					typename impl::state_param state ,
					typename impl::data_param data ) const
			{
//				cout << "Plus transform" << endl;
				return Grammar()( proto::left( expr ) , state , data )
                     + Grammar()( proto::right( expr ) , state , data );
			}
		};
	};


	template< typename Grammar >
	struct minus_transform : proto::transform< minus_transform< Grammar > >
	{
		template< typename Expr , typename State , typename Data >
		struct impl : proto::transform_impl< Expr , State , Data >
		{
			typedef double result_type;

			result_type operator ()(
					typename impl::expr_param expr ,
					typename impl::state_param state ,
					typename impl::data_param data ) const
			{
//				cout << "Minus transform" << endl;
				return Grammar()( proto::left( expr ) , state , data )
                     - Grammar()( proto::right( expr ) , state , data );
			}
		};
	};



	template< typename Grammar >
	struct multiplies_transform : proto::transform< multiplies_transform< Grammar > >
	{
		template< typename Expr , typename State , typename Data >
		struct impl : proto::transform_impl< Expr , State , Data >
		{
			typedef double result_type;

			result_type operator()(
					typename impl::expr_param expr ,
					typename impl::state_param state ,
					typename impl::data_param data ) const
			{
//				cout << "Multiplies transform" << endl;
				typedef typename impl::data data_type;

				double tmp = 0.0;
				for( size_t k=0 ; k<=data.which ; ++k )
				{
					data_type data1( k ,  data.x , data.derivs , data.dt_fac );
					data_type data2( data.which - k ,  data.x , data.derivs , data.dt_fac );

//					tmp += boost::math::binomial_coefficient< double >( data.which , k )
//						* Grammar()( proto::left( expr ) , state , data1 )
//						* Grammar()( proto::right( expr ) , state , data2 );

					tmp += Grammar()( proto::left( expr ) , state , data1 )
						 * Grammar()( proto::right( expr ) , state , data2 );

				}

				return tmp;
			}
		};
	};


	struct terminal_double_transform : proto::transform< terminal_double_transform >
	{
		template< typename Expr , typename State , typename Data >
		struct impl : proto::transform_impl< Expr , State , Data >
		{
			typedef double result_type;

			result_type operator()(
					typename impl::expr_param expr ,
					typename impl::state_param state ,
					typename impl::data_param data ) const
			{
//				cout << "Terminal double transform" << endl;
				return ( data.which == 0 ) ? proto::value( expr ) : 0.0;
			}
		};
	};


	template< typename Index >
	struct terminal_placeholder_transform : proto::transform< terminal_placeholder_transform< Index > >
	{
		template< typename Expr , typename State , typename Data >
		struct impl : proto::transform_impl< Expr , State , Data >
		{
			typedef double result_type;

			result_type operator()(
					typename impl::expr_param expr ,
					typename impl::state_param state ,
					typename impl::data_param data ) const
			{
//				cout << "Terminal placeholder transform" << endl;
				typedef typename impl::expr expr_type;
				typedef typename expr_type::proto_args args_type;
				typedef typename args_type::child0 index_type;
				const size_t index = index_type::value;

				return ( data.which == 0 ) ? data.x[ index ] : data.derivs[ data.which - 1 ][ index ];
			}
		};
	};

	template< typename Grammar >
	struct right_shift_transform : proto::transform< right_shift_transform< Grammar > >
	{
		template< typename Expr , typename State , typename Data >
		struct impl : proto::transform_impl< Expr , State , Data >
		{
			typedef double result_type;

			result_type operator ()(
					typename impl::expr_param expr ,
					typename impl::state_param state ,
					typename impl::data_param data ) const
			{
//				cout << "Right shift transform" << endl;
				return Grammar()( proto::left( expr ) , state , data )
						+ ( data.which == 0 ) ? proto::value( proto::right( expr ) ) : 0.0;
			}
		};
	};

	template< typename Grammar >
	struct left_shift_transform : proto::transform< left_shift_transform< Grammar > >
	{
		template< typename Expr , typename State , typename Data >
		struct impl : proto::transform_impl< Expr , State , Data >
		{
			typedef double result_type;

			result_type operator ()(
					typename impl::expr_param expr ,
					typename impl::state_param state ,
					typename impl::data_param data ) const
			{
//				cout << "Right shift transform" << endl;
				return Grammar()( proto::right( expr ) , state , data )
						+ ( data.which == 0 ) ? proto::value( proto::left( expr ) ) : 0.0;
			}
		};
	};


	template< typename Grammar >
	struct right_scalar_multiplies_transform : proto::transform< right_scalar_multiplies_transform< Grammar > >
	{
		template< typename Expr , typename State , typename Data >
		struct impl : proto::transform_impl< Expr , State , Data >
		{
			typedef double result_type;

			result_type operator ()(
					typename impl::expr_param expr ,
					typename impl::state_param state ,
					typename impl::data_param data ) const
			{
//				cout << "Right scalar multiplies transform" << endl;
				return Grammar()( proto::left( expr ) , state , data ) * proto::value( proto::right( expr ) );
			}
		};
	};

	template< typename Grammar >
	struct left_scalar_multiplies_transform : proto::transform< left_scalar_multiplies_transform< Grammar > >
	{
		template< typename Expr , typename State , typename Data >
		struct impl : proto::transform_impl< Expr , State , Data >
		{
			typedef double result_type;

			result_type operator ()(
					typename impl::expr_param expr ,
					typename impl::state_param state ,
					typename impl::data_param data ) const
			{
//				cout << "Left scalar multiplies transform" << endl;
				return Grammar()( proto::right( expr ) , state , data ) * proto::value( proto::left( expr ) );
			}
		};
	};





	struct taylor_double_terminal : proto::when< proto::terminal< double > , terminal_double_transform > { };

	template< typename Which >
	struct taylor_placeholder_terminal : proto::when< proto::terminal< placeholder< Which > > , terminal_placeholder_transform< Which > > { };

	template< typename Grammar >
	struct taylor_plus : proto::when< proto::plus< Grammar , Grammar > , plus_transform< Grammar > > { };

	template< typename Grammar >
	struct taylor_minus : proto::when< proto::minus< Grammar , Grammar > , minus_transform< Grammar > > { };

	template< typename Grammar >
	struct taylor_multiplies : proto::when< proto::multiplies< Grammar , Grammar > , multiplies_transform< Grammar > >  { };

	template< typename Grammar >
	struct taylor_shift :
		proto::or_<
			proto::when< proto::plus< Grammar , proto::terminal< double > > , right_shift_transform< Grammar > > ,
			proto::when< proto::plus< proto::terminal< double > , Grammar > , left_shift_transform< Grammar > >
		> { };

	template< typename Grammar >
	struct taylor_scalar_multiplies :
		proto::or_<
			proto::when< proto::multiplies< Grammar , proto::terminal< double > > , right_scalar_multiplies_transform< Grammar > > ,
			proto::when< proto::multiplies< proto::terminal< double > , Grammar > , left_scalar_multiplies_transform< Grammar > >
		> { };





	struct taylor_transform :
	proto::or_
	<
		taylor_shift< taylor_transform > ,
		taylor_scalar_multiplies< taylor_transform > ,
		taylor_double_terminal ,
		taylor_placeholder_terminal< proto::_ > ,
		taylor_plus< taylor_transform > ,
		taylor_minus< taylor_transform > ,
		taylor_multiplies< taylor_transform >
	>
	{ };


	template< class System , class State , size_t Order >
	struct eval_derivs
	{
		typedef State state_type;
		typedef taylor_context< state_type , Order > taylor_context_type;
		typedef typename taylor_context_type::deriv_type deriv_type;

		System m_sys;
		taylor_context_type m_data;

		eval_derivs( System sys , const State &x , deriv_type &derivs , double &dt_fac , size_t which )
		: m_sys( sys ) , m_data( which , x , derivs , dt_fac ) { }

		template< class Index >
		void operator()( Index )
		{
			typedef typename fusion::result_of::at< System , Index >::type expr_type;
			const expr_type &expr = boost::fusion::at< Index >( m_sys );

			double deriv = taylor_transform()( expr , 0.0 , m_data );
			m_data.derivs[ m_data.which ][ Index::value ] = m_data.dt_fac / double( m_data.which + 1 ) * deriv;
		}
	};

	template< class System , class State , size_t Order >
	eval_derivs< System , State , Order > make_eval_derivs( System sys , const State &x , std::tr1::array< State , Order > &derivs , double &dt_fac , size_t i )
	{
		return eval_derivs< System , State , Order >( sys , x , derivs , dt_fac , i );
	}
}

template< size_t N , size_t Order >
class taylor
{
public:

	static const size_t dim = N;
	static const size_t order_value = Order;
	typedef double value_type;
	typedef value_type time_type;
	typedef unsigned short order_type;

	typedef std::tr1::array< value_type , dim > state_type;
	typedef state_type deriv_type;
	typedef std::tr1::array< state_type , order_value > derivs_type;

	taylor( value_type rel_error = 1.0e-14 )
	: m_derivs() , m_dt_fac( 1.0 ) , m_rel_error( rel_error ) { }

    order_type order( void ) const
    {
    	return order_value;
    }

    order_type stepper_order( void ) const
    {
    	return order_value;
    }

    order_type error_order( void ) const
    {
    	return order_value - 1;
    }




	template< class System >
	void try_step( System sys , state_type &x , time_type &t , time_type &dt )
	{
		BOOST_STATIC_ASSERT(( boost::fusion::traits::is_sequence< System >::value ));
		BOOST_STATIC_ASSERT(( size_t( boost::fusion::result_of::size< System >::value ) == dim ));

		try_step( sys , x , t , x , dt );
	}

	template< class System >
	void try_step( System sys , const state_type &in , time_type &t , state_type &out , time_type &dt )
	{
		BOOST_STATIC_ASSERT(( boost::fusion::traits::is_sequence< System >::value ));
		BOOST_STATIC_ASSERT(( size_t( boost::fusion::result_of::size< System >::value ) == dim ));

		eval_derivs( sys , in , m_derivs , m_dt_fac );

		double max_error = 0.0;
		for( size_t i=0 ; i<dim ; ++i )
		{
			double error = std::abs( m_derivs[order_value-1][i] ) / ( std::abs( in[i] ) + 1.0e-35 );
			max_error = std::max( error , max_error );
		}

		dt = pow( m_rel_error / max_error , 1.0 / double( order_value ) );
//		clog << dt << tab << max_error << tab << in[0] << tab << in[1] << tab << in[2] << tab;
//		clog << m_derivs[0][0] << tab << m_derivs[0][1] << tab << m_derivs[0][2] << tab << m_dt_fac << endl;

		for( size_t i=0 ; i<dim ; ++i )
		{
			value_type tmp = 0.0;
			for( size_t k=0 ; k<order_value ; ++k )
				tmp = dt * ( tmp + m_derivs[order_value-1-k][i] );
			out[i] = in[i] + tmp;
		}

//		clog << dt << tab << dt0 << tab << max_error << tab << m_dt_fac;
//		for( size_t j=0 ; j<dim ; ++j ) clog << tab << out[j];
//		clog << endl;

//		clog << endl;
//		for( size_t i=0 ; i<4 ; ++i )
//		{
//			for( size_t j=0 ; j<dim ; ++j )
//				clog << m_derivs[i][j] << "\t";
//			clog << endl;
//		}
//		clog << endl;

		dt *= m_dt_fac;
		t += dt;
	}

	const derivs_type& get_last_derivs( void ) const
	{
		return m_derivs;
	}

	template< class System >
	void eval_derivs( System sys , const state_type &in , derivs_type &der , double &dt_fac ) const
	{
		const double min_error = 1.0e-19;
		const double max_error = 1.0e19;
		const double min_fac = 1.5;
		const double max_fac = 0.6;

		for( size_t i=0 ; i<order_value ; ++i )
		{
			boost::mpl::for_each< boost::mpl::range_c< size_t , 0 , dim > >( make_eval_derivs( sys , in , der , dt_fac , i ) );

//			clog << i << tab << "Deriv : ";
//			for( size_t j=0 ; j<dim ; ++j )
//				clog << tab << der[i][j];
//			clog << endl;

			/*
			 * OK # Fehler bestimmen
			 * OK # Fehler grenzen pruefen
			 * OK  # Falls ja
			 *     OK # dt_fac aendern
			 *     OK # schon berechnete Ableitungen skalieren
			 */
			while( true )
			{
				double err = 0.0;
				for( size_t j=0 ; j<dim ; ++j ) err += std::abs( der[i][j] );
//				clog << i;
//				for( size_t j=0 ; j<dim ; ++j ) clog << tab << der[i][j];
//				clog << tab << err << endl;

				if( err < min_error )
				{
//					clog << i << tab << "min_error : " << err << tab << dt_fac << endl;
					scale_derivs( der , i , min_fac );
					dt_fac *= min_fac;
					continue;
				}
				if( err > max_error )
				{
//					clog << i << tab << "max_error : " << err << tab << dt_fac << endl;
					scale_derivs( der , i , max_fac );
					dt_fac *= max_fac;
					continue;
				}
				break;
			}

		}

	}

	void scale_derivs( derivs_type &der , size_t order , double fac ) const
	{
		double scale = 1.0;
		for( size_t i=0 ; i<=order ; ++i )
		{
			scale *= fac;
			for( size_t j=0 ; j<dim ; ++j )
				der[i][j] *= scale;
		}
	}



private:

	derivs_type m_derivs;
	time_type m_dt_fac;
	value_type m_rel_error;

};





} // namespace odeint
} // namespace numeric
} // namespace boost


#endif /* BOOST_NUMERIC_ODEINT_STEPPER_TAYLOR_HPP_ */
