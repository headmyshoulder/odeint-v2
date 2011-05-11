from os import popen
from os import system
from os.path import isfile

toolset = "intel"

#bin_path = "bin/gcc-4.3/release/"
bin_path = "bin/intel-linux/release/"
#bin_path = "bin\\msvc-9.0express\\release\\threading-multi\\"
#extension = ".exe"
extension = ""

bins = [ "odeint_rk4_lorenz" , "odeint_rk4_lorenz_def_alg" , "generic_rk4_lorenz" , "nr_rk4_lorenz" , "gsl_rk4_lorenz" , "rt_generic_rk4_lorenz" ] 
# "odeint_rk54ck" , "odeint_rk54ck_def_alg" , "generic_rk54ck" , "nr_rk54ck" , "gsl_rk54ck" ]

results = []

print "Performance tests for " , bin_path
print

for bin in bins:
	system( "bjam --toolset=" + toolset + " -a " + bin );
	if isfile( bin_path + bin + extension):
		print "Running" , bin
		res = popen( bin_path+bin+extension ).read()
		print bin , res
		results.append( res )
	else:
		results.append( " -- " )

print "Results from" , bin_path
print

for i in range(len(bins)):
	print bins[i] , results[i]
