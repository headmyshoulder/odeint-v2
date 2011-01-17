//============================================================================
// Name        : odeint_ref_problem.cpp
// Author      : Karsten Ahnert
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <ostream>
#include <tr1/array>

#include <boost/function.hpp>
#include <boost/utility.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include <boost/timer.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/detail/iterator_algebra.hpp>
#include <boost/numeric/odeint/container_traits.hpp>


#define tab "\t"

using namespace std;

typedef boost::accumulators::accumulator_set<
    double , boost::accumulators::stats<
		boost::accumulators::tag::mean ,
		boost::accumulators::tag::variance
		>
    > accumulator_type;

typedef boost::timer timer_type;

ostream& operator<<( ostream& out , accumulator_type &acc )
{
    out << boost::accumulators::mean( acc ) << tab;
//    out << boost::accumulators::variance( acc ) << tab;
    return out;
}




class explicit_euler
{
public:

	typedef double value_type;
	typedef std::tr1::array< double , 3 > state_type;
	typedef boost::function3< void , const state_type& , state_type& , value_type > func_type;

	void do_step( func_type system , state_type &x , value_type t , value_type dt )
	{
		system( x , m_dxdt , t );
		for( size_t i=0 ; i<x.size() ; ++i ) x[i] += dt * m_dxdt[i];
	};

private:

	state_type m_dxdt;
};

class explicit_euler_2
{
public:

	typedef double value_type;
	typedef std::tr1::array< double , 3 > state_type;
	typedef boost::numeric::odeint::container_traits< state_type > traits_type;

	template< class System >
	void do_step( System system , state_type &x , value_type t , value_type dt )
	{
		typename boost::unwrap_reference< System >::type sys = system;
		sys( x , m_dxdt , t );
        // x = x + dt*dxdt
        boost::numeric::odeint::detail::it_algebra::increment(
        		traits_type::begin(x) ,
        		traits_type::end(x) ,
        		traits_type::begin(m_dxdt) ,
        		dt );
	};

private:

	state_type m_dxdt;
};


typedef explicit_euler::value_type value_type;
typedef explicit_euler::state_type state_type;

ostream& operator<<( ostream& out , const state_type &x )
{
	out << x[0] << tab << x[1] << tab << x[2];
	return out;
}

const value_type sigma = 10.0;
const value_type R = 28.0;
const value_type b = 8.0 / 3.0;


void lorenz( const state_type &x , state_type &dxdt , value_type t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

struct lorenz_class
{
	void operator()( const state_type &x , state_type &dxdt , value_type t )
	{
	    dxdt[0] = sigma * ( x[1] - x[0] );
	    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
	    dxdt[2] = x[0]*x[1] - b * x[2];
	}
};


int main()
{
	explicit_euler euler1;
	boost::numeric::odeint::stepper_euler< state_type > euler2;
	explicit_euler_2 euler3;
	const value_type dt = 0.01;
	lorenz_class lorenz_c;

	timer_type timer;
	accumulator_type acc1 , acc2 , acc3 , acc4 , acc5 , acc6 , acc7 , acc8 , acc9 , acc10;
	double elapsed;
	const size_t num_of_steps = 1024 * 1024 * 32;

	while( true )
	{
		state_type x = {{ drand48() , drand48() , drand48() }};
		state_type x1( x ) , x2( x ) , x3( x ) , x4( x ) , x5( x ) , x6( x ) , x7( x ) , x8( x ) , x9( x ) , x10( x );

		timer.restart();
		for( size_t i = 0 ; i<num_of_steps ; ++i )
			euler1.do_step( lorenz , x1 , 0.0 , dt );
		elapsed = timer.elapsed();
		acc1( elapsed );
//		clog << "Method 1 with function pointer : " << elapsed << tab << x1 << endl;

		timer.restart();
		for( size_t i = 0 ; i<num_of_steps ; ++i )
			euler1.do_step( lorenz_class() , x2 , 0.0 , dt );
		elapsed = timer.elapsed();
		acc2( elapsed );
//		clog << "Method 1 with class : " << elapsed << tab << x2 << endl;

		timer.restart();
		for( size_t i = 0 ; i<num_of_steps ; ++i )
			euler1.do_step( lorenz_c , x3 , 0.0 , dt );
		elapsed = timer.elapsed();
		acc3( elapsed );
//		clog << "Method 1 with class and copying : " << elapsed << tab << x3 << endl;

		timer.restart();
		for( size_t i = 0 ; i<num_of_steps ; ++i )
			euler1.do_step( boost::ref( lorenz_c ) , x4 , 0.0 , dt );
		elapsed = timer.elapsed();
		acc4( elapsed );
//		clog << "Method 1 with class and boost::ref : " << elapsed << tab << x4 << endl;

		timer.restart();
		for( size_t i = 0 ; i<num_of_steps ; ++i )
			euler2.do_step( lorenz , x5 , 0.0 , dt );
		elapsed = timer.elapsed();
		acc5( elapsed );
//		clog << "Method 2 with function pointer : " << elapsed << tab << x5 << endl;

		timer.restart();
		for( size_t i = 0 ; i<num_of_steps ; ++i )
			euler2.do_step( lorenz_c , x6 , 0.0 , dt );
		elapsed = timer.elapsed();
		acc6( elapsed );
//		clog << "Method 2 with class reference : " << elapsed << tab << x6 << endl;


		timer.restart();
		for( size_t i = 0 ; i<num_of_steps ; ++i )
			euler3.do_step( lorenz , x7 , 0.0 , dt );
		elapsed = timer.elapsed();
		acc7( elapsed );
//		clog << "Method 3 with function pointer : " << elapsed << tab << x6 << endl;

		timer.restart();
		for( size_t i = 0 ; i<num_of_steps ; ++i )
			euler3.do_step( lorenz_c , x8 , 0.0 , dt );
		elapsed = timer.elapsed();
		acc8( elapsed );
//		clog << "Method 3 with class : " << elapsed << tab << x6 << endl;

		timer.restart();
		for( size_t i = 0 ; i<num_of_steps ; ++i )
			euler3.do_step( lorenz_class() , x9 , 0.0 , dt );
		elapsed = timer.elapsed();
		acc9( elapsed );
//		clog << "Method 3 with class : " << elapsed << tab << x6 << endl;

		timer.restart();
		for( size_t i = 0 ; i<num_of_steps ; ++i )
			euler3.do_step( boost::ref( lorenz_c ) , x10 , 0.0 , dt );
		elapsed = timer.elapsed();
		acc10( elapsed );
//		clog << "Method 3 with class reference: " << elapsed << tab << x6 << endl;





		cout << x1 << tab << x2 << tab << x3 << tab << x4 << tab << x5 << tab << x6 << x7 << tab << x8 << tab << x9 << tab << x10 << endl;


		clog << acc1 << tab << acc2 << tab << acc3 << tab << acc4 << tab << acc5 << tab << acc6 << tab << acc7 << tab << acc8 << tab << acc9 << tab << acc10 << endl;
	}

	return 0;
}
