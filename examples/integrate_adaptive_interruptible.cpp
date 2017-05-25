/*
 * integrate_adaptive_intr.cpp
 *
 * This example is a unit test for showing how to run interruptible adaptive ODE45
 * with some points forced to the system. This example forces all points starting 
 * from 0 with 0.1 interval. The output resuls show that timing points such as 
 * 0, 0.1, 0.2, 0.3, 0.4, ... are forced to the solver. This point forcing is useful
 * for some techniques such as Model Predictive Control (MPC)
 *
 * Copyright 2011-2012 Karsten Ahnert
 * Copyright 2011-2013 Mario Mulansky
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#include <iostream>
#include <boost/numeric/odeint.hpp>
// #include "libs/odeint.hpp"

using namespace std;
using namespace boost::numeric::odeint;

typedef double state_type;

class CSystem
{
public:
    
    void operator() ( const state_type &x , state_type &dxdt , const double t)
    {
        dxdt=-x+0.2*input(t);
    }

    double input(const double t) const
    {
        return t>0?1:0;
    }

    void observer( const state_type &x , const double t )
    {
        cout<< t << "\t" << x <<endl;
    }

    double interrupt(const state_type &x, const double t)
    {
        return t+0.1;
    }

    void solve()
    {
        typedef runge_kutta_dopri5<state_type> stepper_type;
        state_type x;
        x= 0.0;
        std::function<void(const state_type &,const double)> my_observer = [&](const state_type &x,const double t){observer(x,t);};
        integrate_adaptive(make_controlled(1E-10,1E-10,stepper_type()),
                            *this,x,0.0,3.0,0.1,my_observer);
    }

    void solve_interruptible()
    {
        typedef runge_kutta_dopri5<state_type> stepper_type;
        state_type x;
        x= 0.0;
        std::function<void(const state_type &,const double)> my_observer = [&](const state_type &x,const double t){observer(x,t);};
        std::function<double(const state_type &,const double)> my_interrupt = [&](const state_type &x,const double t)->double{return interrupt(x,t);};
        integrate_adaptive_interruptible(make_controlled(1E-10,1E-10,stepper_type()),
                            *this,x,0.0,3.0,0.1,my_interrupt,my_observer);
    }

};

int main()
{
    CSystem my_system;
    cout<<"*********** solve normal ***************"<<endl;
    my_system.solve();
    cout<<"*********** solve interruptable ***************"<<endl;
    my_system.solve_interruptible();
    cout<<"*********** end ***************"<<endl;
}

/*
Output results:

*********** solve normal ***************
0   0
3.88978e-07 7.07049e-08
7.77957e-07 1.48501e-07
2.52836e-06 4.98581e-07
1.04052e-05 2.07393e-06
4.58508e-05 9.16287e-06
0.000205356 4.106e-05
0.000923131 0.000184534
0.00415312  0.000828894
0.0186881   0.0037029
0.0702934   0.0135759
0.121899    0.0229524
0.173504    0.0318573
0.22511 0.0403143
0.278719    0.0486496
0.332329    0.0565498
0.385939    0.0640376
0.44152 0.0713884
0.497101    0.0783417
0.552682    0.0849191
0.61034 0.0913668
0.667998    0.0974532
0.725655    0.103199
0.785508    0.108822
0.84536 0.114119
0.905212    0.119109
0.967391    0.123985
1.02957 0.128568
1.09175 0.132874
1.15641 0.137077
1.22107 0.141017
1.28573 0.14471
1.35304 0.148309
1.42035 0.151674
1.48766 0.15482
1.55781 0.157881
1.62797 0.160734
1.69813 0.163395
1.77135 0.165979
1.84457 0.168381
1.91779 0.170614
1.99433 0.172779
2.07087 0.174785
2.1474  0.176643
2.22755 0.178442
2.30769 0.180102
2.38784 0.181634
2.47191 0.183115
2.55599 0.184477
2.64007 0.185729
2.72846 0.186936
2.81685 0.188041
2.90835 0.189087
2.99985 0.190041
3   0.190043
*********** solve interruptible ***************
0   0
3.88978e-07 7.07049e-08
7.77957e-07 1.48501e-07
2.52836e-06 4.98581e-07
1.04052e-05 2.07393e-06
4.58508e-05 9.16287e-06
0.000205356 4.106e-05
0.000923131 0.000184534
0.00415312  0.000828894
0.0186881   0.0037029
0.0702934   0.0135759
0.1 0.0190325
0.152388    0.028269
0.2 0.0362538
0.253359    0.0447621
0.3 0.0518364
0.354614    0.0597112
0.4 0.065936
0.455885    0.0732227
0.5 0.0786939
0.55717 0.0854344
0.6 0.0902377
0.65847 0.0964715
0.7 0.100683
0.759787    0.106447
0.8 0.110134
0.861121    0.115462
0.9 0.118686
0.962475    0.123611
1   0.126424
1.06385 0.130975
1.1 0.133426
1.16524 0.137631
1.2 0.139761
1.26666 0.143646
1.3 0.145494
1.3681  0.149082
1.4 0.150681
1.46957 0.153995
1.5 0.155374
1.57106 0.158435
1.6 0.159621
1.67258 0.162448
1.7 0.163463
1.77413 0.166074
1.8 0.16694
1.87571 0.169351
1.9 0.170086
1.97731 0.172312
2   0.172933
2.07895 0.174988
2.1 0.175509
2.18061 0.177406
2.2 0.177839
2.28231 0.17959
2.3 0.179948
2.37959 0.181482
2.4 0.181856
2.48571 0.183347
2.5 0.183583
2.56432 0.184606
2.6 0.185145
2.68887 0.186408
2.7 0.186559
2.75008 0.187215
2.8 0.187838
2.89217 0.188909
2.9 0.188995
2.93522 0.189376
3   0.190043
*********** end ***************

*/