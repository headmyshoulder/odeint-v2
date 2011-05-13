/*
 * main.cpp
 *
 *  Created on: Apr 2, 2011
 *      Author: karsten
 */

#include <iostream>
#include <fstream>

#include "taylor.hpp"

#include <boost/fusion/include/make_vector.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/assign.hpp>
#include <boost/timer.hpp>

#include <boost/mpl/range_c.hpp>
#include <boost/mpl/for_each.hpp>

#define tab "\t"

template< typename T , size_t N >
std::ostream& operator<<( std::ostream& out , const std::tr1::array< T , N > &x )
{
	if( !x.empty() ) out << x[0];
	for( size_t i=1 ; i<x.size() ; ++i ) out << "\t" << x[i];
	return out;
}

typedef std::tr1::array< double , 3 > state_type;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;
namespace mpl = boost::mpl;

using namespace std;
using namespace boost::numeric::odeint;
using namespace boost::assign;

using boost::numeric::odeint::taylor_adf::arg1;
using boost::numeric::odeint::taylor_adf::arg2;
using boost::numeric::odeint::taylor_adf::arg3;

struct run
{
	vector< double > m_eps_abs_values;
	vector< double > m_eps_rel_values;
	double m_t_end;

	run( const vector< double > &eps_abs_values , const vector< double > &eps_rel_values , double t_end )
	: m_eps_abs_values( eps_abs_values ) , m_eps_rel_values( eps_rel_values ) , m_t_end( t_end ){ }

	template< class Order >
	void operator()( Order ) const
	{
		static const size_t order = Order::value;
		typedef boost::numeric::odeint::taylor< 3 , order > taylor_type;

		boost::timer timer;

		char str[512] = "";
		sprintf( str , "dat/lorenz_taylor_%02d.dat" , int( order ) );
		ofstream fout( str );

		for( size_t i=0 ; i<m_eps_abs_values.size() ; ++i )
		{
			for( size_t j=0 ; j<m_eps_rel_values.size() ; ++j )
			{
				double eps_abs = m_eps_abs_values[i];
				double eps_rel = m_eps_rel_values[j];

				taylor_type taylor( eps_rel , eps_abs );

				state_type x = {{ 10.0 , 10.0 , 10.0 }};

				timer.restart();
				size_t steps_taylor = 0;
				double t = 0.0 , dt = 1.0;
				while( t < m_t_end )
				{
					taylor.try_step(
						fusion::make_vector
						(
							sigma * ( arg2 - arg1 ) ,
							R * arg1 - arg2 - arg1 * arg3 ,
							arg1 * arg2 - b * arg3
						) , x , t , dt );
					steps_taylor++;
				}
				double time_taylor = timer.elapsed();

				fout << i << tab << j << tab << eps_abs << tab << eps_rel << tab;
				fout << steps_taylor << tab << time_taylor;
				fout << endl;
			}
		}
	}
};




int main( int argc , char **argv )
{
	if( argc != 2 )
	{
		cerr << "usage taylor_lorenz_eval_adaptive_statistics t_end" << endl;
		return -1;
	}
	double t_end = atof( argv[1] );

	vector< double > eps_abs_values;
	vector< double > eps_rel_values;

	eps_abs_values += 1.0e1 , 1.0 , 1.0e-1 , 1.0e-2 , 1.0e-3 , 1.0e-4 , 1.0e-5 , 1.0e-6 , 1.0e-7 , 1.0e-8 , 1.0e-9 , 1.0e-10 , 1.0e-11 , 1.0e-12 , 1.0e-13 , 1.0e-14;
	eps_rel_values += 1.0e-1 , 1.0e-2 , 1.0e-3 , 1.0e-4 , 1.0e-5 , 1.0e-6 , 1.0e-7 , 1.0e-8 , 1.0e-9 , 1.0e-10 , 1.0e-11 , 1.0e-12 , 1.0e-13 , 1.0e-14;


	typedef mpl::range_c< size_t , 5 , 30 > order_values;
	mpl::for_each< order_values >( run( eps_abs_values , eps_rel_values , t_end ) );



	return 0;
}
