
/*
 [auto_generated]
 boost/numeric/odeint/iterator/const_step_ode_iterator.hpp

 [begin_description]
 Iterator for iterating throught the solution of an ODE with constant step size.
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_ITERATOR_CONST_STEP_ODE_ITERATOR_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_ITERATOR_CONST_STEP_ODE_ITERATOR_HPP_INCLUDED

# include <boost/iterator/iterator_facade.hpp>


namespace boost {
namespace numeric {
namespace odeint {


    template< class Stepper , class System , class StepperTag = typename Stepper::stepper_category > 
    class const_step_iterator;



    // Specilization for steppers and error steppers
    template< class Stepper , class System >
    class const_step_iterator< Stepper , System , stepper_tag > : public boost::iterator_facade
    <
        const_step_iterator< Stepper , System , stepper_tag > ,
        typename Stepper::state_type const ,
        boost::single_pass_traversal_tag
    >
    {
    private:

        typedef Stepper stepper_type;
        typedef System system_type;
        typedef typename stepper_type::state_type state_type;
        typedef typename stepper_type::time_type time_type;

    public:
   
        const_step_iterator( stepper_type stepper , system_type sys , state_type &s , time_type t , time_type dt )
            : m_stepper( stepper ) , m_system( sys ) , m_state( s ) , m_t( t ) , m_dt( dt ) {}

    private:

        friend class boost::iterator_core_access;

        void increment()
        {
            m_stepper.do_step( m_system , m_state , m_t , m_dt );
            m_t += m_dt;
        }

        bool equal( const_step_iterator const& other ) const
        {
            return m_t > other.m_t;
        }

        const state_type& dereference() const
        {
            return m_state;
        }

        stepper_type m_stepper;
        system_type m_system;
        state_type &m_state;
        time_type m_t;
        time_type m_dt;

    };



    // Specilization for dense output stepper
    template< class Stepper , class System >
    class const_step_iterator< Stepper , System , dense_output_stepper_tag > : public boost::iterator_facade
    <
        const_step_iterator< Stepper , System , dense_output_stepper_tag > ,
        typename Stepper::state_type const ,
        boost::single_pass_traversal_tag
    >
    {
    private:

        typedef Stepper stepper_type;
        typedef System system_type;
        typedef typename stepper_type::state_type state_type;
        typedef typename stepper_type::time_type time_type;

    public:
   
        const_step_iterator( stepper_type stepper , system_type sys , state_type &s , time_type t , time_type dt )
            : m_stepper( stepper ) , m_system( sys ) , m_state( s ) , m_t( t ) , m_dt( dt )
        {
            m_stepper.initialize( m_state , m_t , m_dt );
        }

    private:

        friend class boost::iterator_core_access;

        void increment( void )
        {
            m_t += m_dt;
            while(  m_stepper.current_time() < m_state.second )
                m_stepper.do_step( m_system );
            m_stepper.calc_state( m_t , m_state );
        }

        bool equal( const_step_iterator const& other ) const
        {
            return m_t > other.m_t;
        }

        const state_type& dereference( void ) const
        {
            return m_state;
        }

        stepper_type m_stepper;
        system_type m_system;
        state_type &m_state;
        time_type m_t;
        time_type m_dt;
    };





    /*
     * ToDo : create specializations for controlled steppers and dense output steppers
     */


    template< class Stepper , class System >
    const_step_iterator< Stepper , System > make_const_step_iterator(
        Stepper stepper ,
        System system , 
        typename Stepper::state_type &x ,
        typename Stepper::time_type t ,
        typename Stepper::time_type dt )
    {
        return const_step_iterator< Stepper , System >( stepper , system , x , t , dt );
    }


    template< class Stepper , class System >
    std::pair< const_step_iterator< Stepper , System > , const_step_iterator< Stepper , System > >
    make_const_step_range(
        Stepper stepper ,
        System system , 
        typename Stepper::state_type &x ,
        typename Stepper::time_type t_start ,
        typename Stepper::time_type t_end ,
        typename Stepper::time_type dt )
    {
        return std::make_pair(
            const_step_iterator< Stepper , System >( stepper , system , x , t_start , dt ) ,
            const_step_iterator< Stepper , System >( stepper , system , x , t_end , dt ) );
    }












} // namespace odeint
} // namespace numeric
} // namespace boost



#endif // BOOST_NUMERIC_ODEINT_ITERATOR_CONST_STEP_ODE_ITERATOR_HPP_INCLUDED
