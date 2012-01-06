/*
 * main.cpp
 *
 *  Created on: Apr 2, 2011
 *      Author: karsten
 */

#include <iostream>
#include <fstream>

#include <boost/numeric/odeint/config.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <boost/numeric/odeint/stepper/generation.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>


#include <boost/fusion/include/make_vector.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/assign.hpp>
#include <boost/timer.hpp>

#include <boost/mpl/range_c.hpp>
#include <boost/mpl/for_each.hpp>

#include "taylor.hpp"

#define tab "\t"

template< typename T , size_t N >
std::ostream& operator<<( std::ostream& out , const boost::array< T , N > &x )
{
	if( !x.empty() ) out << x[0];
	for( size_t i=1 ; i<x.size() ; ++i ) out << "\t" << x[i];
	return out;
}

typedef boost::array< double , 3 > state_type;

typedef boost::numeric::odeint::runge_kutta_cash_karp54< state_type > rk54_type;
typedef boost::numeric::odeint::runge_kutta_dopri5< state_type > dopri5_type;
typedef boost::numeric::odeint::runge_kutta_fehlberg78< state_type > fehlberg_type;
typedef boost::numeric::odeint::bulirsch_stoer< state_type > bulirsch_stoer_type;


const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

struct print {
    ofstream &m_out;
    print( ofstream &out )
        : m_out( out )
    { }

    void operator()( const state_type &x , const double t )
    {
        m_out << x[0] << tab << x[1] << tab << x[2] << endl;
    }
};

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


	boost::timer timer;
	ofstream fout1( "dat/lorenz_rk54.dat" );
	ofstream fout2( "dat/lorenz_dopri5.dat" );
	ofstream fout3( "dat/lorenz_fehlberg.dat" );
	ofstream fout4( "dat/lorenz_bs.dat" );

	ofstream wf_out( "dat/bs.dat" );

	for( size_t i=0 ; i<eps_abs_values.size() ; ++i )
	{
		for( size_t j=0 ; j<eps_rel_values.size() ; ++j )
		{
			double eps_abs = eps_abs_values[i];
			double eps_rel = eps_rel_values[j];

			rk54_type rk54;
			dopri5_type dopri5;
			fehlberg_type fehlberg;
			bulirsch_stoer_type bulirsch_stoer( eps_abs , eps_rel );


			state_type x1 = {{ 10.0 , 10.0 , 10.0 }};
			state_type x2 = {{ 10.0 , 10.0 , 10.0 }};
			state_type x3 = {{ 10.0 , 10.0 , 10.0 }};
	        state_type x4 = {{ 10.0 , 10.0 , 10.0 }};

			timer.restart();
			size_t steps_rk54 = integrate_adaptive( make_controlled( eps_abs , eps_rel , rk54 ) , lorenz , x1 , 0.0 , t_end , 0.1 );
			double time_rk54 = timer.elapsed();

			fout1 << i << tab << j << tab << eps_abs << tab << eps_rel << tab;
			fout1 << steps_rk54 << tab << time_rk54 << tab;
			fout1 << endl;



			timer.restart();
			size_t steps_dopri5 = integrate_adaptive( make_controlled( eps_abs , eps_rel , dopri5 ) , lorenz , x2 , 0.0 , t_end , 0.1 );
			double time_dopri5 = timer.elapsed();

			fout2 << i << tab << j << tab << eps_abs << tab << eps_rel << tab;
			fout2 << steps_dopri5 << tab << time_dopri5 << tab;
			fout2 << endl;



	        timer.restart();
	        size_t steps_fehlberg = integrate_adaptive( make_controlled( eps_abs , eps_rel , fehlberg ) , lorenz , x3 , 0.0 , t_end , 0.1 );
	        double time_fehlberg = timer.elapsed();

	        fout3 << i << tab << j << tab << eps_abs << tab << eps_rel << tab;
	        fout3 << steps_fehlberg << tab << time_fehlberg << tab;
	        fout3 << endl;



	        timer.restart();
	        size_t steps_bs = integrate_adaptive( bulirsch_stoer , lorenz , x4 , 0.0 , t_end , 0.1 ); //, print(wf_out) );
	        double time_bs = timer.elapsed();

	        fout4 << i << tab << j << tab << eps_abs << tab << eps_rel << tab;
	        fout4 << steps_bs << tab << time_bs << tab;
	        fout4 << endl;

            clog << "done" << endl;
		}
		fout1 << endl;
		fout2 << endl;
		fout3 << endl;
		fout4 << endl;
	}




	return 0;
}
