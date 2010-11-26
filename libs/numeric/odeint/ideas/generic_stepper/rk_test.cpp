#include "runge_kutta_stepper.hpp"

int main( int argc , char **argv )
{
	typedef runge_kutta_stepper< 4 > rk_type;
	typedef rk_type::coef_a_type coef_a_type;
	typedef rk_type::coef_b_type coef_b_type;
	typedef rk_type::coef_c_type coef_c_type;

	const boost::array< double , 1 > a1 = {{ 5.1 }};
	const boost::array< double , 2 > a2 = {{ 5.2 , 5.3 }};
	const boost::array< double , 3 > a3 = {{ 5.4 , 5.5 , 5.6 }};

	const coef_a_type a = fusion::make_vector( a1 , a2 , a3 );
	const coef_b_type b = {{ 6.1 , 6.2 , 6.3 , 6.4 }};
	const coef_c_type c = {{ 1.1 , 1.2 , 1.3 , 1.4 }};

	rk_type rk( a , b , c );
	rk.print_vals();

}
