#! /usr/bin/python

import fnmatch
import os

preamble_template = "/*\n\
 FILENAME\n\
 \n\
 [begin_description]\n\
 \n\
 [end_description]\n\
 \n\
 Copyright 2009-2011 Karsten Ahnert\n\
 Copyright 2009-2011 Mario Mulansky\n\
 \n\
 Distributed under the Boost Software License, Version 1.0.\n\
 (See accompanying file LICENSE_1_0.txt or\n\
 copy at http://www.boost.org/LICENSE_1_0.txt)\n\
*/\n\
\n\
"

def header_guard( file_desc ):
	return file_desc[2].upper().replace( "/" , "_" ).replace( "." , "_" ) + "_INCLUDED"

def insert_preamble( content , file_desc ):
	preamble = preamble_template.replace( "FILENAME" , file_desc[2] )
	ret = preamble + "\n"
	ret += "#ifndef " + header_guard( file_desc ) + "\n"
	ret += "#define " + header_guard( file_desc ) + "\n\n"
	ret += content
	ret += "#endif // " + header_guard( file_desc ) + "\n"
	return ret

	
def glob_dir( dir , pattern , exclude = " " ):
	matches = []
	for root , dirnames, filenames in os.walk( dir ):
		if root.find( exclude ) == -1:
			for filename in fnmatch.filter( filenames , pattern ):
				matches.append( [ root , filename , os.path.join(root, filename) ] )
	return matches

files = glob_dir( "boost" , "*.hpp" , "stepper" )

for f in files:
	# print f + " " + header_guard( f )
	file = open( f[2] )
	content = file.read()
	content = insert_preamble( content , f )
	file.close()
	
	file = open( f[2] , "w" )
	file.write( content )
	# file.close()