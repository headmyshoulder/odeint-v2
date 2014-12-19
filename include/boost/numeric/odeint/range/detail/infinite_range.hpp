/*
  [auto_generated]
  boost/numeric/odeint/range/detail/infinite_range.hpp

  [begin_description]
  tba.
  [end_description]

  Copyright 2009-2012 Karsten Ahnert
  Copyright 2009-2012 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_RANGE_DETAIL_INFINITE_RANGE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_RANGE_DETAIL_INFINITE_RANGE_HPP_INCLUDED


namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< typename Range , typename Value >
class infinite_range
{
public:

    using range_type = Range;
    using value_type = Value;
    using reference = value_type const&;
    using base_type = infinite_range< range_type , value_type >;


    struct iterator : public std::iterator< std::input_iterator_tag , value_type >
    {
        iterator( base_type *rng ) : m_rng( rng ) {}

        iterator operator++( void )
        {
            m_rng->increment();
            return *this;
        }

        iterator operator++( int )
        {
            iterator orig = *this;
            ++(*this);
            return orig;
        }

        reference operator*( void ) const
        {
            return m_rng->dereference();
        }

        bool operator==( iterator const& other ) const
        {
            return m_rng == other.m_rng;
        }

        bool operator!=( iterator const& other ) const
        {
            return m_rng != other.m_rng;
        }

        base_type * m_rng;
    };

    iterator begin( void )
    {
        return iterator( this );
    }

    iterator end( void )
    {
        return iterator( nullptr );
    }

protected:

    range_type const& range( void ) const
    {
        return *static_cast< range_type const* >( this );
    }

    range_type& range( void )
    {
        return *static_cast< range_type* >( this );
    }

    reference dereference( void ) const
    {
        return range().dereference();
    }

    void increment( void )
    {
        range().increment();
    }
};



} // namespace detail
} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_RANGE_DETAIL_INFINITE_RANGE_HPP_INCLUDED
