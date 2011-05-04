from os import popen
from os.path import isfile

bin_path = "bin/gcc-4.6/release/"
#bin_path = "bin/intel-linux/release/"

#bins = [ "odeint_rk4" , "odeint_rk4_def_alg" , "generic_rk4" , "nr_rk4" , "gsl_rk4" , "rt_generic_rk4" , 
bins = [ "odeint_rk54ck" , "odeint_rk54ck_def_alg" , "generic_rk54ck" , "nr_rk54ck" , "gsl_rk54ck" ]
results = []

print "Performance tests for " , bin_path
print

for bin in bins:
	if isfile( bin_path + bin ):
		print "Running" , bin
		res = popen( bin_path+bin ).read()
		print bin , res
		results.append( res )
	else:
		results.append( " -- " )

print "Results from" , bin_path
print

for i in range(len(bins)):
	print bins[i] , results[i]
