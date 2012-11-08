/*
 [auto_generated]
 boost/numeric/odeint/stepper/base/algebra_stepper_base.hpp

 [begin_description]
 Base class for all steppers with an algebra and operations.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_BASE_ALGEBRA_STEPPER_BASE_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_BASE_ALGEBRA_STEPPER_BASE_HPP_INCLUDED


namespace boost {
namespace numeric {
namespace odeint {

/**
 * \class algebra_stepper_base
 * \brief Base class for all steppers with algebra and operations.
 *
 * This class serves a base class for all steppers with algebra and operations. It holds the
 * algebra and provides access to the algebra.  The operations are not instantiated, since they are 
 * static classes inside the operations class.
 *
 * \tparam Algebra The type of the algebra. Must fullfil the Algebra Concept, at least partially to work
 * with the stepper.
 * \tparam Operations The type of the operations. Must fullfil the Operations Concept, at least partially 
 * to work with the stepper.
 */
template< class Algebra , class Operations >
class algebra_stepper_base
{
public:

    typedef Algebra algebra_type;
    typedef Operations operations_type;

    /**
     * \brief Constructs a algebra_stepper_base and creates the algebra. This constructor can be used as a default
     * constructor if the algebra has a default constructor.
     * \param algebra The algebra_stepper_base stores and uses a copy of algebra.
     */
    algebra_stepper_base( const algebra_type &algebra = algebra_type() )
    : m_algebra( algebra ) { }

    /**
     * \return A reference to the algebra which is held by this class.
     */
    algebra_type& algebra()
    {
        return m_algebra;
    }

    /**
     * \return A const reference to the algebra which is held by this class.
     */
    const algebra_type& algebra() const
    {
        return m_algebra;
    }

protected:

    algebra_type m_algebra;
};

} // odeint
} // numeric
} // boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_BASE_ALGEBRA_STEPPER_BASE_HPP_INCLUDED
