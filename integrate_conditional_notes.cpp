/*

= no own stepper =

no own stepper, since this stepper cannot be used in the integrate routines at least not without throwing an exception!

*/

integrate_const( cond_stepper , sys , x , 0.0 , 10.0 , 0.1 ); // can not stop at t = 5.0;


/*

= idea for conditional integrate =

question: will stepper perform const step size integration?
*/

integrate_conditional( stepper , sys , x , 0.0 , controller );


void integrate_conditional( S1 stepper , S2 sys , S3 x , T1 t , T2 dt , C controller )
{
    controller.init( stepper , sys , x , t , dt );
    while( !controller.stop( x , t ) )
    {
        controller.do_step( stepper , sys , x , t , dt );
    }
    controller.exit( sys , x , t , dt );
}

template< class Time >
struct adaptive_controller
{
    typedef Time time_type;

    template< class Stepper , class Sys , class State >
    void init( Stepper &stepper , Sys sys , const State &s , time_type t , time_type dt )
    {
        typedef typename Stepper::stepper_category stepper_category;
        init_impl( stepper , s , t , dt , stepper_category() );
    }

    template< class State >
    bool stop( const State &s , time_type &t )
    {
        return t > m_t_end;
    }

    template< class Stepper , class Sys , class State >
    void do_step( Stepper &stepper , Sys sys , State &x , time_type t , time_type dt )
    {
        typedef typename Stepper::stepper_category stepper_category;
        do_step_impl( stepper , sys , x , t , dt , stepper_category() )
    }

    template< class Stepper , class Sys , class State >
    void exit( Stepper &stepper , Sys sys , State &x , time_type t , time_type dt )
    {
        typedef typename Stepper::stepper_category stepper_category;
        exit_impl( stepper , x , t );
    }

    





    template< class Stepper , class Sys , class State >
    void do_step_impl( Stepper &stepper , Sys sys , State &x , time_type &t , time_type &dt , stepper_tag )
    {
        stepper.do_step( sys , x , t , dt );
        t += dt;
    }

    template< class Stepper , class Sys , class State >
    void do_step_impl( Stepper &stepper , Sys sys , State &x , time_type &t , time_type &dt , controlled_stepper_tag )
    {
        size_t trials = 0;
        controlled_step_result res = success;
        do
        {
            res = stepper.try_step( sys , x , t , dt );
            ++trials;
        }
        while( ( res == fail ) && ( trials < max_attempts ) );
        if( trials == max_attempts ) throw std::overflow_error( "adaptive_controller::do_step_impl : too much iterations" );

    }

    template< class Stepper , class Sys , class State >
    void do_step_impl( Stepper &stepper , Sys sys , State &x , time_type &t , time_type &dt , dense_output_stepper_tag )
    {
        stepper.do_step( sys );
        t = stepper.current_time();
    }



    bool m_stop;
    Time m_t_end;
};



/*

 * implement integrate_conditional
 * implement adaptive_controller
 * implement const_step_controller
 * implement adaptive_touch_controller
 * implement const step touch controller
 * unit testing
 * examples
 * documentation
 */
