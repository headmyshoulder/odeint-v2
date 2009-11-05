#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

#include <boost/lambda/lambda.hpp>

#define tab "\t" 

using namespace std;
using namespace boost::lambda;

template< class Iterator1 , class Iterator2 >
void increment( Iterator1 first1 , Iterator1 last1 , Iterator2 first2 )
{
    while( first1 < last1 )
	(*first1++) += (*first2++);
}

template< class Iterator1 , class Iterator2 , class T >
void increment( Iterator1 first1 , Iterator1 last1 , Iterator2 first2 , T val )
{
    while( first1 < last1 )
	(*first1++) += val * (*first2++);
}

template< class Iterator1 , class Iterator2 , class Operation >
void increment_op( Iterator1 first1 , Iterator1 last1 , Iterator2 first2 , Operation op )
{
    while( first1 < last1 )
	(*first1++) += op(*first2++);
}

template< class Iterator1 , class Iterator2 >
void decrement( Iterator1 first1 , Iterator1 last1 , Iterator2 first2 )
{
    while( first1 < last1 )
	(*first1++) -= (*first2++);
}

template< class Iterator1 , class Iterator2 , class T >
void decrement( Iterator1 first1 , Iterator1 last1 , Iterator2 first2 , T val )
{
    while( first1 < last1 )
	(*first1++) -= val * (*first2++);
}




const size_t n = 100000;
const size_t cycles = 500000;

int main( int argc , char **argv )
{
    vector<double> eins( n ) , zwei( n );
    generate( eins.begin() , eins.end() , drand48 );
    generate( zwei.begin() , zwei.end() , drand48 );
    vector<double> eins1 , zwei1;

    cout << "Berechne a += b, mit " << endl;
    cout << "1.  a += b" << endl;
    cout << "2.  a += 1.0 * b" << endl;
    cout << "3.  a -= (-1.0)*b" << endl;
    cout << "4.  a += lambda( _1 )" << endl;
    cout << "5.  a += lambda( 1.0 *_1)" << endl;


    clock_t st1 , et1 , st2 , et2 , st3 , et3 , st4 , et4 , st5 , et5;


    eins1 = eins;
    zwei1 = zwei;
    st1 = clock();
    for( size_t i=0 ; i<cycles ; ++i )
	increment( eins1.begin() , eins1.end() , zwei1.begin() );
    et1 = clock();
    double time1 = double( et1 - st1 ) / double( CLOCKS_PER_SEC );
    cout << "1." << tab << time1 << endl;

    eins1 = eins;
    zwei1 = zwei;
    st2 = clock();
    for( size_t i=0 ; i<cycles ; ++i )
	increment( eins1.begin() , eins1.end() , zwei1.begin() , 1.0 );
    et2 = clock();
    double time2 = double( et2 - st2 ) / double( CLOCKS_PER_SEC );
    cout << "2." << tab << time2 << endl;


    eins1 = eins;
    zwei1 = zwei;
    st3 = clock();
    for( size_t i=0 ; i<cycles ; ++i )
	decrement( eins1.begin() , eins1.end() , zwei1.begin() , -1.0 );
    et3 = clock();
    double time3 = double( et3 - st3 ) / double( CLOCKS_PER_SEC );
    cout << "3." << tab << time3 << endl;

    eins1 = eins;
    zwei1 = zwei;
    st4 = clock();
    for( size_t i=0 ; i<cycles ; ++i )
	increment_op( eins1.begin() , eins1.end() , zwei1.begin() , _1 );
    et4 = clock();
    double time4 = double( et4 - st4 ) / double( CLOCKS_PER_SEC );
    cout << "4." << tab << time4 << endl;

    eins1 = eins;
    zwei1 = zwei;
    st5 = clock();
    for( size_t i=0 ; i<cycles ; ++i )
	increment_op( eins1.begin() , eins1.end() , zwei1.begin() , 1.0*_1 );
    et5 = clock();
    double time5 = double( et5 - st5 ) / double( CLOCKS_PER_SEC );
    cout << "5." << tab << time5 << endl;


    cout << "1." << tab << "2." << tab << "3." << tab << "4." << tab << "5." << endl;
    cout << time1 << tab << time2 << tab << time3 << tab << time4 << tab << time5 << endl;

    return 0;
}


/*
int main( int argc , char **argv )
{
    const size_t n = 1000000;
    vector<double> eins( n ) , zwei( n ) , drei( n );
    generate( eins.begin() , eins.end() , drand48 );
    generate( zwei.begin() , zwei.end() , drand48 );
    generate( drei.begin() , drei.end() , drand48 );

    increment( eins.begin() , eins.end() , zwei.begin() , 0.5*_1 );
    increment( eins.begin() , eins.end() , zwei.begin() , drei.begin() , _1 + 0.5*_2 );

    return 0;
}
*/
