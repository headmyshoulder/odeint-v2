from os import popen
from os import system
from os.path import isfile
from pylab import *

#toolset = "gcc-4.4"
#toolset = "intel"
#toolset = "msvc"
toolset = "msvc-10.0"

#bin_path = "bin/gcc-4.4/release/"
#bin_path = "bin/intel-linux/release/"
bin_path = "bin\\msvc-10.0\\release\\" #threading-multi\\"
extension = ".exe"
#extension = ""

bins = [ "odeint_rk4_lorenz_array" , "odeint_rk4_lorenz_range" , "generic_odeint_rk4_lorenz" , "nr_rk4_lorenz" ,  "rt_generic_rk4_lorenz" ]

results = []

print "Performance tests for " , bin_path
print

for bin in bins:
	#system( "bjam toolset=" + toolset + " -a " + bin );
	if isfile( bin_path + bin + extension):
		print "Running" , bin
		res = popen( bin_path+bin+extension ).read()
		print bin , res
		results.append( res )
	else:
		print "no executable found:" , bin_path + bin + extension
		results.append( 0 )

print "Results from" , bin_path
print

for i in range(len(bins)):
	print bins[i] , results[i]

res = zeros( len(results) )
base = float(results[0])
for i in xrange(len(results)):
	res[i] = base/float(results[i])

bar_width = 0.6

figure(1)
title("Runge-Kutta 4" , fontsize=20)
bar( arange(5) , res , bar_width , color='blue' , linewidth=4 , edgecolor='blue' , ecolor='red') #, elinewidth=2, ecolor='red' )
xlim( -0.5 , 4.5+bar_width )
xticks( arange(5)+bar_width/2 , ('array' , 'range' , 'generic' , 'NR' , 'rt gen' ) )
ylabel('Performance in %' , fontsize=20)

show()