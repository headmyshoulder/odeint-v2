#include <iostream>
#include <cmath>

#include <boost/timer.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

#include <boost/numeric/odeint/stepper/error_checker_explicit.hpp>
#include <boost/numeric/odeint/stepper/generic_controlled_stepper.hpp>
#include <boost/numeric/odeint/stepper/generic_controlled_stepper_explicit.hpp>
#include <boost/numeric/odeint/stepper/controller/default_controller.hpp>
#include <boost/numeric/odeint/stepper/controller/pi_controller.hpp>


#define tab "\t"

namespace boost {
namespace numeric {
namespace odeint {




/*
 * move this class into controller_runge_kutta.hpp
 */
template< class Stepper >
class controlled_runge_kutta_explicit :
    public generic_controlled_stepper<
        Stepper ,
        error_checker_explicit_new< typename Stepper::value_type , typename Stepper::algebra_type , typename Stepper::operations_type > ,
        default_controller ,
        typename Stepper::resizer_type ,
        typename Stepper::stepper_category >
{
public:

    typedef generic_controlled_stepper<
        Stepper ,
        error_checker_explicit_new< typename Stepper::value_type , typename Stepper::algebra_type , typename Stepper::operations_type > ,
        default_controller ,
        typename Stepper::resizer_type ,
        typename Stepper::stepper_category > base_type;
    typedef error_checker_explicit_new< typename Stepper::value_type , typename Stepper::algebra_type , typename Stepper::operations_type > error_checker_type;
    typedef typename Stepper::value_type value_type;

    controlled_runge_kutta_explicit( const Stepper &stepper = Stepper() ,
            const value_type eps_abs = static_cast< value_type >( 1.0e-6 ) ,
            const value_type eps_rel = static_cast< value_type >( 1.0e-6 ) ,
            const value_type a_x = static_cast< value_type >( 1.0 ) ,
            const value_type a_dxdt = static_cast< value_type >( 1.0 ) )
    : base_type( stepper ,
                 error_checker_type( base_type::stepper().algebra() , eps_abs , eps_rel , a_x , a_dxdt )
               ) { }
};


}
}
}




typedef boost::array< double , 3 > state_type;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;

    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = -b * x[2] + x[0] * x[1];
}


int main( int argc , char **argv )
{
    using namespace std;
    using namespace boost::numeric::odeint;

    generic_controlled_stepper<
        runge_kutta_cash_karp54< state_type > ,
        error_checker_explicit< double , range_algebra , default_operations > ,
        default_controller ,
        initially_resizer ,
        explicit_error_stepper_tag > stepper1;

    controlled_runge_kutta< runge_kutta_cash_karp54< state_type > > stepper2;
    controlled_runge_kutta_explicit< runge_kutta_cash_karp54< state_type > > stepper3;

    generic_controlled_stepper<
        runge_kutta_cash_karp54< state_type > ,
        error_checker_explicit< double , range_algebra , default_operations > ,
        pi_controller< double > ,
        initially_resizer ,
        explicit_error_stepper_tag > stepper4;


    state_type x1 = {{ 10.0 , 10.0 , 10.0 }};
    state_type x2 = x1;
    state_type x3 = x1;
    state_type x4 = x1;

    boost::timer timer;
    const double t_max = 100000.0;

    timer.restart();
    size_t steps1 = integrate_adaptive( stepper1 , lorenz , x1 , 0.0 , t_max , 0.1 );
    double t1 = timer.elapsed();
    cout << steps1 << tab << t1 << tab << x1[0] << tab << x1[1] << tab << x1[2] << endl;

    timer.restart();
    size_t steps2 = integrate_adaptive( stepper2 , lorenz , x2 , 0.0 , t_max , 0.1 );
    double t2 = timer.elapsed();
    cout << steps2 << tab << t2 << tab << x2[0] << tab << x2[1] << tab << x2[2] << endl;

    timer.restart();
    size_t steps3 = integrate_adaptive( stepper3 , lorenz , x3 , 0.0 , t_max , 0.1 );
    double t3 = timer.elapsed();
    cout << steps3 << tab << t3 << tab << x3[0] << tab << x3[1] << tab << x3[2] << endl;

    timer.restart();
    size_t steps4 = integrate_adaptive( stepper4 , lorenz , x4 , 0.0 , t_max , 0.1 );
    double t4 = timer.elapsed();
    cout << steps4 << tab << t4 << tab << x4[0] << tab << x4[1] << tab << x4[2] << endl;




    return 0;
}
