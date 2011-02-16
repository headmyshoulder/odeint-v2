#! /usr/bin/python

import os
import glob
import subprocess


topleveldir = "../../../"

filename = []
filename += glob.glob( topleveldir + "boost/numeric/odeint/*.hpp" )
filename += glob.glob( topleveldir + "boost/numeric/odeint/algebra/*.hpp" )
filename += glob.glob( topleveldir + "boost/numeric/odeint/algebra/detail/*.hpp" )
filename += glob.glob( topleveldir + "boost/numeric/odeint/algebra/external/*.hpp" )
filename += glob.glob( topleveldir + "boost/numeric/odeint/stepper/*.hpp" )
filename += glob.glob( topleveldir + "boost/numeric/odeint/stepper/base/*.hpp" )
filename += glob.glob( topleveldir + "boost/numeric/odeint/stepper/detail/*.hpp" )


test_prog = '#include <INCLUDE>\n\
\n\
int main( void )\n\
{\n\
    return 0;\n\
}\n\
'

count = 0
for fn in filename :
    fn = fn[len(topleveldir) : ]
    prog = test_prog.replace( "INCLUDE" , fn )

    file = open( "test.cc" , "w" )
    file.write( prog )
    file.close

    cmd = ["g++" , "-I"+topleveldir , "-I$BOOST_ROOT" , "-c" , "test.cc" ]
    
    print "Test " + fn + " : "
#    print prog
    print "Commandline : " + subprocess.list2cmdline( cmd )
    
    p = subprocess.Popen( subprocess.list2cmdline( cmd ) , shell=True)
    sts = os.waitpid(p.pid, 0)[1]

    print sts
    print ""

