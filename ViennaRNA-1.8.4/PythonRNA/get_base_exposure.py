#!/usr/bin/python

import RNA
from base_exposure import get_base_exposure
from sys import argv

sequence = argv[ 1 ]
temperature = 37 #273.15 + 37
if len( argv ) > 2:  temperature = float( argv[2] )

base_exposure = get_base_exposure( sequence, temperature )
for i in base_exposure: print i,
print
