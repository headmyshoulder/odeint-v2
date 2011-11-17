/*
 * mtl_bindings.hpp
 *
 *  Created on: Oct 26, 2011
 *      Author: mario
 */

#ifndef MTL_BINDINGS_HPP_
#define MTL_BINDINGS_HPP_


//[ mtl_bindings
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/odeint/util/state_wrapper.hpp>

namespace boost { namespace numeric { namespace odeint {

template< class Value , class Parameters >
struct is_resizeable< mtl::dense_vector< Value , Parameters > >
{ // declare resizeablility
    typedef boost::true_type type;
    const static bool value = type::value;
};

template< class Value , class Parameters > //, class Value2 , class Parameters2 >
struct same_size_impl< mtl::dense_vector< Value , Parameters > , mtl::dense_vector< Value , Parameters > >
{ // define how to check size
    static bool same_size( const mtl::dense_vector< Value , Parameters > &v1 ,
                           const mtl::dense_vector< Value , Parameters > &v2 )
    {
        return mtl::size( v1 ) == mtl::size( v2 );
    }
};

template< class Value , class Parameters >
struct resize_impl< mtl::dense_vector< Value , Parameters > , mtl::dense_vector< Value , Parameters > >
{ // define how to resize
    static void resize( mtl::dense_vector< Value , Parameters > &v1 ,
                        const mtl::dense_vector< Value , Parameters > &v2 )
    {
        v1.change_dim( mtl::size( v2 ) );
    }
};

} } }

//]

#endif /* MTL_BINDINGS_HPP_ */
