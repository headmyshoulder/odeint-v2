#include <algorithm>
#include <vector>
#include <iostream>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include <boost/timer.hpp>

#define tab "\t" 

using namespace std;
using namespace boost::accumulators;

using boost::make_zip_iterator;
using boost::make_tuple;
using boost::tuples::get;



typedef accumulator_set<
    double , stats< tag::mean , tag::variance >
    > accumulator_type;

typedef boost::timer timer_type;

ostream& operator<<( ostream& out , accumulator_type &acc )
{
    out << boost::accumulators::mean( acc ) << tab;
//    out << boost::accumulators::variance( acc ) << tab;
    return out;
}



// computes y = x1 + alpha2*x2 + alpha3*x3
template < class Iterator1 , class Iterator2 , class Iterator3 , class T >
inline void scale_sum( Iterator1 first1 , Iterator1 last1 ,
                       T a1 , Iterator2 first2 ,
                       T a2 , Iterator3 first3 )
{
    while( first1 != last1 )
        (*first1++) = a1 * (*first2++) + a2 * (*first3++);
}

template < class Iterator1 , class Iterator2 , class Iterator3 , class Operator >
void transform3_1( Iterator1 first1 , Iterator1 last1,
		   Iterator2 first2 , Iterator3 first3 ,
		   Operator op )
{
    while( first1 != last1 ) op( *first1++ , *first2++ , *first3++ );
}

template < class Iterator1 , class Iterator2 , class Iterator3 , class Operator >
void transform3_2( Iterator1 first1 , Iterator1 last1,
                 Iterator2 first2 , Iterator3 first3 ,
                 Operator op )
{
    for( ; first1 != last1 ; ++first1 , ++first2 , ++first3 )
	op( *first1 , *first2 , *first3 );
}

template < class Iterator1 , class Iterator2 , class Iterator3 , class Operator >
void transform3_3( Iterator1 first1 , Iterator1 last1,
                 Iterator2 first2 , Iterator3 first3 ,
                 Operator op )
{
    while( first1 != last1 )
    {
	op( *first1 , *first2 , *first3 );
	++first1;
	++first2;
	++first3;
    }
}

struct scale_sum_op
{
    double a1 , a2;

    scale_sum_op( double _a1 , double _a2 )
        : a1( _a1 ) , a2( _a2 ) {}

    template< class T >
    void operator()( T &x1 , T &x2 , T &x3 )
    {
        x1 = a1 * x2 + a2 * x3;
    }
};

struct scale_sum_op_tuple
{
    double a1 , a2;
    scale_sum_op_tuple( double _a1 , double _a2 )
        : a1( _a1 ) , a2( _a2 ) {}

    template< class Tuple >
    void operator()( const Tuple &t ) const
    {
	get<0>(t) = a1 * get<1>(t) + a2 * get<2>(t);
    }
};
                 

int main( int argc , char **argv )
{
    const size_t n = 1024 * 1024;
    const size_t num_of_iterations = 16;
    double a1 = 0.25 , a2 = 1.25;


    accumulator_type acc1 , acc2 , acc3 , acc4 , acc5;
    timer_type timer;


    size_t count = 0;
    clog.precision(4);
    while( true )
    {
        vector< double > org1( n ) , org2( n ) , org3( n );
        generate( org1.begin() , org1.end() , drand48 );
        generate( org2.begin() , org2.end() , drand48 );
        generate( org3.begin() , org3.end() , drand48 );

        vector< double > v1( org1 ) , v2( org2 ) , v3( org3 );
        vector< double > w1( org1 ) , w2( org2 ) , w3( org3 );
        vector< double > u1( org1 ) , u2( org2 ) , u3( org3 );
        vector< double > x1( org1 ) , x2( org2 ) , x3( org3 );

        timer.restart();
        for( size_t i=0 ; i<num_of_iterations ; ++i )
            scale_sum( v1.begin() , v1.end() , a1 , v2.begin() , a2 , v3.begin() );
        double res1 = accumulate( v1.begin() , v1.end() , 0.0 );
        acc1( timer.elapsed() );

        timer.restart();
        for( size_t i=0 ; i<num_of_iterations ; ++i )
            transform3_1( w1.begin() , w1.end() , w2.begin() , w3.begin() ,
                        scale_sum_op( a1 , a2 ) );
        double res2 = accumulate( w1.begin() , w1.end() , 0.0 );
        acc2( timer.elapsed() );


	timer.restart();
	for( size_t i=0 ; i<num_of_iterations ; ++i )
	    std::for_each(
		make_zip_iterator( make_tuple( u1.begin() , u2.begin() , u3.begin() ) ) ,
		make_zip_iterator( make_tuple( u1.end() , u2.end() , u3.end() ) ) ,
		scale_sum_op_tuple( a1 , a2 )
		);
	double res3 = accumulate( u1.begin() , u1.end() , 0.0 );
	acc3( timer.elapsed() );

	timer.restart();
	for( size_t i=0 ; i<num_of_iterations ; ++i )
	    transform3_2( x1.begin() , x1.end() , x2.begin() , x3.begin() ,
			  scale_sum_op( a1 , a2 ) );
	double res4 = accumulate( x1.begin() , x1.end() , 0.0 );
	acc4( timer.elapsed() );

	timer.restart();
	for( size_t i=0 ; i<num_of_iterations ; ++i )
	    transform3_3( x1.begin() , x1.end() , x2.begin() , x3.begin() ,
			  scale_sum_op( a1 , a2 ) );
	double res5 = accumulate( x1.begin() , x1.end() , 0.0 );
	acc5( timer.elapsed() );

	if( ( res1 != res2 ) ||
	    ( res2 != res3 ) ||
	    ( res3 != res4 ) ||
	    ( res4 != res5 ) )
	    clog << "error" << endl;

        ++count;

        clog << count << tab;
        clog << acc1 << tab;
        clog << acc2 << tab;
	clog << acc3 << tab;
	clog << acc4 << tab;
	clog << acc5 << tab;
/*	clog << tab;
        clog << res1 << tab;
        clog << res2 << tab;
	clog << res3 << tab;
	clog << res4 << tab;
	clog << res5 << tab;*/

        clog << endl;
    }

    return 0;
}






/*
  Results, average time in

  scale_sum transform3_1 for_each(zip_iterator) transform3_2 transform3_3

  kink with gcc-4.2

  0.153           0.1521          0.151           0.1483          0.1472

  photon with gcc-4.4, kann verfaelscht sein, weil viel laeuft

  0.163           0.1606          0.1646          0.1696          0.1694

  quark with gcc-4.4, kann verfaelscht sein

  0.1574          0.1629          0.1686          0.1731          0.1732

  at home with gcc-4.4

  0.1004          0.1001          0.1021          0.1009          0.1009 
*/
  
