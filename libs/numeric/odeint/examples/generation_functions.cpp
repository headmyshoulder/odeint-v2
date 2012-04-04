/*
 libs/numeric/odeint/examples/stochastic_euler.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Stochastic euler stepper example and Ornstein-Uhlenbeck process

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <array>

#include <boost/numeric/odeint.hpp>

typedef std::array< double , 1 > state_type;

using namespace boost::numeric::odeint;

int main( int argc , char **argv )
{
    {
        typedef runge_kutta_dopri5< state_type > stepper_type;

        //[ generation_functions_syntax_auto
        auto stepper1 = make_controlled( 1.0e-6 , 1.0e-6 , stepper_type() );
        auto stepper2 = make_dense_output( 1.0e-6 , 1.0e-6 , stepper_type() );
        //]

        //[ generation_functions_syntax_result_of
        result_of::make_controlled< stepper_type >::type stepper3 = make_controlled( 1.0e-6 , 1.0e-6 , stepper_type() );
        result_of::make_dense_output< stepper_type >::type stepper4 = make_dense_output( 1.0e-6 , 1.0e-6 , stepper_type() );
        //]
    }
    return 0;
}
