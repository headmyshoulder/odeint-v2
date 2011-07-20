#include <vector>

using namespace std;

struct rt_algebra
{
	template< typename T , size_t dim >
	inline static void foreach( boost::array< T , dim > & x_tmp , const boost::array< T , dim > &x ,
				const vector< double > &a ,
				const boost::array< T , dim > *k_vector , const double dt )
	{
		for( size_t i=0 ; i<dim ; ++i )
        {
            x_tmp[i] = x[i];
            for( size_t j = 0 ; j<a.size() ; ++j )
                x_tmp[i] += a[j]*dt*k_vector[j][i];
        }
	}
};