/*
*/

#include <cmath>

#include <boost/array.hpp>

template< size_t N >
struct phase_lattice
{
	typedef double value_type;
	typedef boost::array< value_type , N > state_type;

	value_type m_epsilon; 
	state_type m_omega;

	phase_lattice() : m_epsilon( 6.0/(N*N) ) // should be < 8/N^2 to see phase locking
	{
		for( size_t i=1 ; i<N-1 ; ++i )
			m_omega[i] = m_epsilon*(N-i);
	}

	void inline operator()( const state_type &x , state_type &dxdt , const double t ) const
	{
		dxdt[0] = m_omega[0] + sin( x[1] - x[0] );

		for( size_t i=1 ; i<N-1 ; ++i )
			dxdt[i] = m_omega[i] + sin( x[i] - x[i-1] ) + sin( x[i+1] - x[i] );

		dxdt[N-1] = m_omega[N-1] + sin( x[N-1] - x[N-2] );
	}

};
