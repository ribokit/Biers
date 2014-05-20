#!/usr/bin/python

from base_exposure import *
from sys import argv

import sys
sys.path.append('/Users/rhiju/python/')
from parse_options import parse_options

offset = parse_options( argv, 'offset', -1 )
mutpos = parse_options( argv, 'mutpos', [-1] )
temperature = parse_options( argv, 'temperature', 37.0 );
only_library1 = parse_options( argv, 'only_library1', 0 );

#print 'OFFSET', offset

sequence = argv[1]
[tapestry,mutations] = get_tapestry( sequence, mutpos, offset, temperature, only_library1 )

for i in range( len( mutations ) ):
    print mutations[i],
    for m in range( len(tapestry[i]) ):
        print tapestry[i][m],
    print
