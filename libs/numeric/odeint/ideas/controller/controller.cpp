/*
 * TODO:
 *
 * * Put in appropriate header files
 *
 * * implement controllers
 *   * default
 *   * pi control
 *
 * * implement error checkers
 *   * different norms
 *
 * * implement overloads for
 *   * explicit fsal
 *   * implicit (stepper)
 *
 * * implement try_step overloads for
 *   * different variants of try_step
 *   * solving the forwarding problem
 *
 * * get_value / set_value (default_operations) nach util auslagern
 *
 * * testing
 */

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

#include "error_checker_explicit.hpp"
#include "generic_controlled_stepper.hpp"
#include "generic_controlled_stepper_explicit.hpp"

#define tab "\t"

namespace boost {
namespace numeric {
namespace odeint {

class default_controller
{
public:

    template< class Value , class Time , class Order >
    controlled_step_result operator()( Value error , Time &t , Time &dt , Order stepper_order , Order error_order )
    {
        using std::max;
        using std::min;
        using std::pow;

        if( error > 1.0 )
        {
            // error too large - decrease dt ,limit scaling factor to 0.2 and reset state
            dt *= max( 0.9 * pow( error , -1.0 / ( Value( error_order ) - 1.0 ) ) , 0.2 );
            return fail;
        }
        else
        {
            if( error < 0.5 )
            {
                //error too small - increase dt and keep the evolution and limit scaling factor to 5.0
                t += dt;
                dt *= min( 0.9 * pow( error , -1.0 / Value( stepper_order ) ) , 5.0 );
                return success;
            }
            else
            {
                t += dt;
                return success;
            }
        }

    }
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

    state_type x1 = {{ 10.0 , 10.0 , 10.0 }};
    state_type x2 = x1;

    boost::timer timer;

    timer.restart();
    size_t steps1 = integrate_adaptive( stepper1 , lorenz , x1 , 0.0 , 10000000.0 , 0.1 );
    double t1 = timer.elapsed();

    timer.restart();
    size_t steps2 = integrate_adaptive( stepper2 , lorenz , x2 , 0.0 , 10000000.0 , 0.1 );
    double t2 = timer.elapsed();

    cout << steps1 << tab << t1 << tab << x1[0] << tab << x1[1] << tab << x1[2] << endl;
    cout << steps2 << tab << t2 << tab << x2[0] << tab << x2[1] << tab << x2[2] << endl;

    return 0;
}
