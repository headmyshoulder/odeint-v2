#! /usr/bin/python

import os
import glob

initial_comment = '/* FILENAME header file\n\
\n\
 Copyright 2010 Karsten Ahnert\n\
 Copyright 2010 Mario Mulansky\n\
\n\
 Distributed under the Boost Software License, Version 1.0.\n\
 (See accompanying file LICENSE_1_0.txt or\n\
 copy at http://www.boost.org/LICENSE_1_0.txt)\n\
*/\n\
\n\
\n\
\n\
'
print initial_comment


topleveldir = "../../../"

filenames = []
filenames += glob.glob( topleveldir + "boost/numeric/odeint/*.hpp" )
filenames += glob.glob( topleveldir + "boost/numeric/odeint/algebra/*.hpp" )
filenames += glob.glob( topleveldir + "boost/numeric/odeint/algebra/detail/*.hpp" )
filenames += glob.glob( topleveldir + "boost/numeric/odeint/algebra/external/*.hpp" )
filenames += glob.glob( topleveldir + "boost/numeric/odeint/stepper/*.hpp" )
filenames += glob.glob( topleveldir + "boost/numeric/odeint/stepper/base/*.hpp" )
filenames += glob.glob( topleveldir + "boost/numeric/odeint/stepper/detail/*.hpp" )


to_find = [ "#include" , "namespace" , "#ifndef" ]

def find_first_of( str , to_find ) :
    indices = []
    for sub in to_find :
        indices.append( str.find( sub ) )
    while( indices.count( -1 ) > 0 ):
        indices.remove( -1 )
    if( len( indices ) == 0 ):
        return -1
    min_index = min( indices )
    return min_index

for fn in filenames :
    content = open( fn ).read()
    start = find_first_of( content , to_find )

    if( start == -1 ) :
        print "Nothing to replace found in " + fn
        continue
    
    if( start == 0 ):
        print "In " + fn + " nothing to replace was found!"
    else :
        print "In " + fn + " the following string was found to be replaced : "
        print content[0:start]
    answer = raw_input( "Should this part be replaced with the Copyright comment? " )
    print answer , answer.lower()
    if( ( answer.lower() == "yes" ) or ( answer.lower() == "y" ) ):
        print content.replace( content[0:start] , initial_comment.replace( "FILENAME" , fn[len(topleveldir):] ) )
    print ""


 
