/*
 boost header: numeric/odeint/hamiltonian_stepper_euler.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_HAMILTONIAN_STEPPER_EULER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_HAMILTONIAN_STEPPER_EULER_HPP_INCLUDED

#include <stdexcept>

#include <boost/numeric/odeint/detail/iterator_algebra.hpp>
#include <boost/numeric/odeint/container_traits.hpp>

namespace boost {
namespace numeric {
namespace odeint {

    template<
        class Container ,
        class Time = double ,
        class Traits = container_traits< Container >
        >
    class hamiltonian_stepper_euler
    {
        // provide basic typedefs
    public:

        typedef unsigned short order_type;
        typedef Time time_type;
        typedef Traits traits_type;
        typedef typename traits_type::container_type container_type;
        typedef typename traits_type::value_type value_type;
        typedef typename traits_type::iterator iterator;
        typedef typename traits_type::const_iterator const_iterator;

    private:

	container_type m_dpdt;



    public:

	template< class CoordinateFunction >
	void do_step( CoordinateFunction &qfunc ,
		      container_type &q ,
		      container_type &p ,
		      time_type dt )
        {
	    if( !traits_type::same_size( q , p ) )
	    {
		std::string msg( "hamiltonian_stepper_euler::do_step(): " );
		msg += "q and p have different sizes";
		throw std::invalid_argument( msg );
	    }

            traits_type::adjust_size( p , m_dpdt );

	    detail::it_algebra::increment( traits_type::begin(q) ,
					   traits_type::end(q) ,
					   traits_type::begin(p) ,
					   dt );
	    qfunc( q , m_dpdt );
	    detail::it_algebra::increment( traits_type::begin(p) ,
					   traits_type::end(p) ,
					   traits_type::begin(m_dpdt) ,
					   dt );
	}

    };


} // namespace odeint
} // namespace numeric
} // namespace boost





#endif //BOOST_NUMERIC_ODEINT_HAMILTONIAN_STEPPER_EULER_HPP_INCLUDED
