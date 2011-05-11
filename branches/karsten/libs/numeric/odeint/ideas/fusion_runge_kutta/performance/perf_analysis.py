from pylab import *

rcParams.update( { 'xtick.labelsize': 20 , 'ytick.labelsize': 20 })


results_rk4 = array( 
[[ 0.82 , 0.59 , 0.60 , 0.60 , 0.54 , 0.54 , 0.81 , 0.98 , 1.14 ] , #odeint
 [ 0.64 , 0.56 , 0.54 , 0.60 , 0.77 , 0.65 , 0.85 , 0.98 , 1.09 ] , #generic
 [ 0.96 , 0.64 , 0.62 , 0.60 , 0.80 , 0.47 , 1.05 , 0.63 , 0.76 ] , #nr
 [ 1.52 , 1.61 , 1.35 , 1.77 , 1.07 , 1.06 , 1.07 , 1.28 , 1.62 ] , #gsl
 [ 1.85 , 1.85 , 1.47 , 2.08 , 1.08 , 1.08 , 1.57 , 1.70 , 2.10 ]]) #rt gen

results_rk54ck = array(
[[ 1.34 , 1.14 , 1.38 , 1.00 , 0.95 , 1.53 , 1.91 ] , #odeint
 [ 1.26 , 1.13 , 1.43 , 1.28 , 1.02 , 1.63 , 1.84 ] , #generic
 [ 2.04 , 1.32 , 1.34 , 1.48 , 1.10 , 1.82 , 1.15 ] , #nr
 [ 2.80 , 2.34 , 2.79 , 1.80 , 1.94 , 1.99 , 2.16 ]]) #gsl

means_rk4 = 100*ones( 5 )
error_rk4 = zeros( 5 )

for i in arange(1,5):
	tmp = results_rk4[0] / results_rk4[i]
	means_rk4[i] = 100*mean( tmp )
	error_rk4[i] = 100*sqrt(var( tmp ))

means_rk54ck = 100*ones( 4 )
error_rk54ck = zeros( 4 )

for i in arange(1,4):
	tmp = results_rk54ck[0] / results_rk54ck[i]
	means_rk54ck[i] = 100*mean( tmp )
	error_rk54ck[i] = 100*sqrt(var( tmp ))

bar_width = 0.6
  
figure(1)
title("Runge-Kutta 4" , fontsize=20)
bar( arange(5) , means_rk4 , bar_width , color='blue' , linewidth=4 , edgecolor='blue' , yerr = error_rk4 , ecolor='red') #, elinewidth=2, ecolor='red' )
xlim( -0.5 , 4.5+bar_width )
xticks( arange(5)+bar_width/2 , ('odeint' , 'generic' , 'NR' , 'GSL' , 'rt gen' ) )
ylabel('Performance in %' , fontsize=20)

figure(2)
title("Runge-Kutta 5(4) Cash-Karp" , fontsize=20)
bar( arange(4) , means_rk54ck , bar_width , color='blue' , linewidth=4 , edgecolor='blue' , yerr = error_rk54ck , ecolor='red' )
xlim( -0.5 , 3.5+bar_width )
xticks( arange(4)+bar_width/2 , ('odeint' , 'generic' , 'NR' , 'GSL' ) )
ylabel('Performance in %' , fontsize=20)

show()
