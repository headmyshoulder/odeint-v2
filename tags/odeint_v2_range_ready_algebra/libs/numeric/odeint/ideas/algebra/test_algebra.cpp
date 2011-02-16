/*
 * test_algebra.cpp
 *
 *  Created on: Feb 12, 2011
 *      Author: karsten
 */

#include <vector>
#include <iostream>

#include "algebra.hpp"

using namespace std;

struct print
{
	template< class T >
	void operator()( T t )
	{
		cout << t << endl;
	}
};

int main( int argc , char **argv )
{
	algebra2 al;
	vector< double > vec( 3 );
	vec[0] = 0.1;
	vec[1] = 0.95;
	vec[2] = 1.23;

	al.for_each1( vec , print() );
	cout << al.tmp << endl;

	algebra1 al2;

//	al2.for_each1()( vec , print() );
//	al.for_each2()( vec , print() );

	return 0;
}
