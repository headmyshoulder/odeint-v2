Iterator versions:

class XYZ_iterator
{
    state_type m_x;
    stepper_type m_stepper;
    system_type m_system;    
};

template< class Stepper , class System >
class XYZ_iterator
{
    typedef Stepper unwrapped_stepper_type;
    typedef System unwrapped_system_type;

    typedef typename odeint::unwrap_reference< unwrapped_stepper_type >::type stepper_type;
    typedef typename stepper_type::state_type state_type;
    typedef typename stepper_type::time_type time_type;
    typedef typename odeint::unwrap_reference< unwrapped_system_type >::type system_type

    public:

    XYZ_iterator( unwrapped_stepper_type stepper , unwrapped_system_type system , state_type &x ,
                  time_type t_start , time_type t_end , time_end dt )
        : m_stepper( stepper ) , m_system( system ) , m_x( &x ) ,
        m_t( t_start ) , m_t_end( t_end ) , m_dt( dt ) , m_is_last( false )
    { }

    state_type *m_x;
    unwrapped_stepper_type m_stepper;
    unwrapped_system_type m_system;
    time_type m_t;
    time_type m_t_end;
    time_type m_dt;
    bool m_is_last;
};



TEST:

test copy constructor
test assignment operator


