import fnmatch
import os


def glob_dir( dir , pattern , exclude = [ ] ):
    matches = []
    for root , dirnames, filenames in os.walk( dir ):
        exclude_dir = False
        for e in exclude:
            if root.find( e ) != -1:
                exclude_dir = True
        if exclude_dir == False: 
            for filename in fnmatch.filter( filenames , pattern ):
                matches.append( [ root , filename , os.path.join(root, filename) ] )
    return matches
