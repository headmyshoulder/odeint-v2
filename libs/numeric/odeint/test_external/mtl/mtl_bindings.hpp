/*
 * mtl_bindings.hpp
 *
 *  Created on: Oct 26, 2011
 *      Author: mario
 */

#ifndef MTL_BINDINGS_HPP_
#define MTL_BINDINGS_HPP_

#include <boost/numeric/mtl/mtl.hpp>

namespace boost { namespace numeric { namespace odeint {

template< class Value1 , class Parameters1 , class Value2 , class Parameters2 >
bool same_size( mtl::dense_vector< Value1 , Parameters1 > &v1 ,
                const mtl::dense_vector< Value2 , Parameters2 > &v2 )
{
    return mtl::size( v1 ) == mtl::size( v2 );
}

template< class Value1 , class Parameters1 , class Value2 , class Parameters2 >
void resize( mtl::dense_vector< Value1 , Parameters1 > &v1 ,
             const mtl::dense_vector< Value2 , Parameters2 > &v2 )
{
    v1.change_dim( mtl::size( v2 ) );
}

} } }

#endif /* MTL_BINDINGS_HPP_ */
