/*
 boost header: numeric/odeint/mtl4_dense2d_container_traits.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Container traits for mtl4::dense2D matrices

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_CONTAINER_TRAITS_MTL4_DENSE2D_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_CONTAINER_TRAITS_MTL4_DENSE2D_HPP_INCLUDED

#include <boost/numeric/odeint/container_traits.hpp>
#include <boost/numeric/mtl/mtl.hpp>

namespace boost {
namespace numeric {
namespace odeint {

    template <typename Value, typename Parameters >
    struct container_traits< mtl::dense2D< Value , Parameters > >
    {

        typedef mtl::matrix::dense2D< Value , Parameters > container_type;

        typedef typename container_type::value_type value_type;

        typedef typename mtl::traits::range_generator<
            mtl::tag::all, container_type >::type cursor;
        typedef typename mtl::traits::value< container_type >::type value_map;
        typedef typename mtl::traits::const_value<
            container_type >::type const_value_map;
        typedef mtl::utilities::iterator_adaptor<
            value_map, cursor, value_type > iterator;
        typedef mtl::utilities::iterator_adaptor<
            const_value_map, cursor, value_type > const_iterator;
        
        static void resize( const container_type &x , container_type &dxdt )
        {
            dxdt = container_type( x.num_rows(), x.num_cols() );
        }
        
        static bool same_size(
            const container_type &x1 ,
            const container_type &x2
            )
        {
            return ( (x1.num_rows() == x2.num_rows()) && (x1.num_cols() == x2.num_cols()));
        }

        static void adjust_size(
                const container_type &x1 ,
                container_type &x2
            )
        {
            if( !same_size( x1 , x2 ) ) resize( x1 , x2 );
        }

        static iterator begin( container_type &x )
        {
            cursor c = mtl::begin< mtl::tag::all >(x);
            value_map v(x);
            return iterator(v, c);
        }

        static const_iterator begin( const container_type &x )
        {
            cursor c = mtl::begin< mtl::tag::all >(x);
            const_value_map v(x);
            return const_iterator(v, c);
        }

        static iterator end( container_type &x )
        {
            cursor c = mtl::end< mtl::tag::all >(x);
            value_map v(x);
            return iterator(v, c);
        }

        static const_iterator end( const container_type &x )
        {
            cursor c = mtl::end< mtl::tag::all >(x);
            const_value_map v(x);
            return const_iterator(v, c);
        }

    };



} // namespace odeint
} // namespace numeric
} // namespace boost



/* Template Specialization to provide += operator for iterator return type */
namespace mtl { namespace utilities { namespace detail {

template <typename Cursor, typename ValueType>
struct iterator_proxy<typename mtl::traits::value< mtl::matrix::dense2D< ValueType > >::type, Cursor, ValueType>
{
    typedef iterator_proxy                    self;
    typedef typename mtl::traits::value< mtl::matrix::dense2D< ValueType > >::type property_map;

    iterator_proxy(property_map map, Cursor cursor) 
        : map(map), cursor(cursor) {}

    operator ValueType() const
    {
        return static_cast<ValueType>(map(*cursor));
    }

    self& operator=(ValueType const& value)
    {
        map(*cursor, value);
        return *this;
    }

    self& operator+=(ValueType const& value)
    {
        map( *cursor, value + static_cast<ValueType>map(*cursor)) );
        return *this;
    }

    property_map           map;
    Cursor                 cursor;
};

}}} // namespace mtl::utilities::detail


// define iterator traits for dense2D iterators
namespace std {

template < typename Cursor , typename Value , typename Parameters>
struct iterator_traits< 
    mtl::utilities::iterator_adaptor< mtl::detail::direct_value< mtl::dense2D< Value , Parameters > >, Cursor, Value >
>
{
    typedef Value value_type;
    typedef int difference_type; // ?
    typedef Value* pointer;
    typedef Value& reference;
    typedef std::random_access_iterator_tag iterator_category;
};


template < typename Cursor , typename Value , typename Parameters>
struct iterator_traits< 
    mtl::utilities::iterator_adaptor< mtl::detail::direct_const_value< mtl::dense2D< Value , Parameters > >, Cursor, Value >
>
{
    typedef Value value_type;
    typedef int difference_type; // ?
    typedef const Value* pointer;
    typedef const Value& reference;
    typedef std::random_access_iterator_tag iterator_category;
};


}


#endif //BOOST_NUMERIC_ODEINT_CONTAINER_TRAITS_MTL4_DENSE2D_HPP_INCLUDED
