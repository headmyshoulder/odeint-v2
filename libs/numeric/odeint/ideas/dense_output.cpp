/*
 * dense_output.cpp
 *
 *  Created on: Sep 29, 2010
 *      Author: karsten
 */

class dopri5_dense_output
{
public:

	template< class System >
	std::pair< time_type , time_type > do_step( System &system );

	state_type get_state( double t );
	void calc_state( double t , state_type &x );
	state_type& get_current_state( void );
};


void test_dense_output
{
	dopri_dense_output dopri;

	dopri.initialize( x0 , t0 , dt0 );

	while( t < tend )
	{
		std::pair< double , double > time_range = dopri.do_step( system );
		state_type x = dopri.get_state( time_range.first + 0.1 );
		dopri.calc_state( time_range.first + 0.1 , x );
	}
}
