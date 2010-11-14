/*
 * main.cpp
 *
 *  Created on: Oct 16, 2010
 *      Author: karsten
 */

#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <typeinfo>

#include <tr1/array>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/inserter.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/push_back.hpp>

#include <boost/fusion/container.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/include/mpl.hpp>
//#include <boost/fusion/adapted/mpl.hpp>
//#include <boost/fusion/container/vector/convert.hpp>
//#include <boost/fusion/include/as_vector.hpp>


#define tab "\t"

using namespace std;
namespace fusion = boost::fusion;
namespace mpl = boost::mpl;

//template< class T , class Constant >
//struct array_wrapper : std::tr1::array< T , Constant::value >
//{
//};

template< class T , class Constant >
struct array_wrapper
{
	typedef typename std::tr1::array< T , Constant::value > type;
};


struct print_xml
{
    template <typename T>
    void operator()(T const& x) const
    {
        std::cout
            << '<' << typeid(x).name() << '>'
//            << x
            << "</" << typeid(x).name() << '>'
            << endl;
    }
};


template< size_t N >
class test_class
{
public:

	typedef typename fusion::result_of::as_vector
	<
		typename mpl::copy
		<
			mpl::range_c< size_t , 1 , N > ,
			mpl::inserter
			<
				mpl::vector0<> ,
				mpl::push_back< mpl::_1 , array_wrapper< double , mpl::_2 > >
			>
		>::type
	>::type coef_a_type;

	typedef std::tr1::array< double , N > coef_c_type;

	const coef_a_type coef_a;
	const coef_c_type coef_c;

	test_class( const coef_a_type &a , const coef_c_type &c ) : coef_a( a ) , coef_c( c ) { }
};



int main( int argc , char **argv )
{
	typedef test_class< 3 > class_type;
	const std::tr1::array< double , 1 > a_1 = {{ 0.5 }};
	const std::tr1::array< double , 2 > a_2 = {{ 2.2 , 3.3 }};
	const class_type::coef_a_type coef_a = fusion::make_vector( a_1 , a_2 );

	test_class< 3 >::coef_a_type coef_a_new;
	fusion::for_each( coef_a_new , print_xml() );


	const class_type::coef_c_type coef_c = {{ 1.0 , 2.0 }};

	class_type my_class( coef_a , coef_c );

	cout << my_class.coef_c[0] << endl;
	cout << my_class.coef_c[1] << endl;

	cout << fusion::at_c< 0 >( my_class.coef_a )[0] << endl;
	cout << fusion::at_c< 1 >( my_class.coef_a )[0] << endl;
	cout << fusion::at_c< 1 >( my_class.coef_a )[1] << endl;

	return 0;
}
