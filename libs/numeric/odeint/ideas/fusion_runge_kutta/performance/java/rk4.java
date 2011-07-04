//package org.apache.commons.math.ode;
//package org.apache.commons.math.ode.nonstiff;

import org.apache.commons.math.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math.ode.FirstOrderIntegrator;
import org.apache.commons.math.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math.ode.DerivativeException;
import org.apache.commons.math.ode.IntegratorException;

public class rk4 {


    private static class Lorenz implements FirstOrderDifferentialEquations 
    {
        private static double sigma = 10.0;
        private static double R = 28.0;
        private static double b = 8.0/3.0;

        public int getDimension() 
        {
            return 3;
        }

        public void computeDerivatives( double t , double[] y , double[] yDot )
        {
            yDot[0] = sigma * ( y[1] - y[0] );
            yDot[1] = R * y[0] - y[1] - y[0] * y[2];
            yDot[2] = y[0]*y[1] - b * y[2];
        }
    }


    public static void main( String[] args )
    {
        final int loops = 20;
        final double dt = 1E-10;
        final int steps = 20000000;
        double mean = 0.0;

        FirstOrderIntegrator rk4 = new ClassicalRungeKuttaIntegrator( dt );
        final FirstOrderDifferentialEquations ode = new Lorenz();

        for( int n=0 ; n<loops+3 ; n++ )
        {
            double[] y = new double[] { 0.0, 1.0 , 1.0 }; // initial state
            long start = System.currentTimeMillis();
            try {
                rk4.integrate(ode, 0.0, y, dt*steps , y);
            } catch(DerivativeException de) {
                System.out.println("wrong exception caught");
            } catch(IntegratorException ie) {
            }
            System.out.println( y[0] + "  " + y[1] + "  " + y[2] );
            long end = System.currentTimeMillis();
            System.out.println( "Time: " + (end-start) + " ms" );
            if( n > 2 )
                mean += end - start;
        }
        System.out.println("Average run time: " + (mean/loops) + " ms" );
    }

}