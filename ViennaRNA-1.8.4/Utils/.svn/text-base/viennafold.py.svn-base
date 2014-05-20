#!/usr/bin/python
import math
import RNA
import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-x', '--bonusmatrix', dest='xfilename', default='',\
                  help='A file that contains the matrix of experimental bonuses, in csv format',\
                  metavar='XMATRIXFILE')
parser.add_option('-s', '--sh', dest='shfilename', default='',\
                  help='A file that contains the shape bonuses, given in RNAstructure format.\
                        These constraints specificially use SHAPE pseudoenergy constraints',\
                  metavar='SHAPEFILE')
parser.add_option('-b', '--si', dest='si', default=-0.8, type='float',\
                  help='Specify an intercept used with SHAPE constraints.\
                        Default is -0.8 kcal/mol.')
parser.add_option('-m', '--sm', dest='sm', default=2.6, type='float',\
                  help='Specify a slope used with SHAPE constraints.\
                        Default is 2.6 kcal/mol.')
(options, args) = parser.parse_args()
#I have to adjust to kcal/mol since Vienna does everything in cal/mol
adjust = 100
def parse_sequence(filename):
    linen = 0
    for l in open(filename).readlines():
	l = l.strip()
	if l[0] != ';':
	    if linen == 0:
		linen += 1
	    else:
		return l.replace(' ','').upper().strip('1')

def mysplit(s, delim):
    res = []
    for ch in s.split(delim):
	if ch:
	    res.append(ch)
    return res

sequence = parse_sequence(args[0])
outfile = open(args[1], 'w')
RNA.cvar.noLonelyPairs = 1
RNA.cvar.pf_scale = 1.0
RNA.fold(sequence)
print 'Folding %s ...' % sequence
if options.xfilename:
    bonusfilename ='/tmp/viennaRNA_bonuses.txt'
    bonusfile = open(bonusfilename, 'w')
    bonuses = []
    xbonuses = [0]*len(sequence)
    for i in range(len(sequence)):
	bonuses.append([0]*len(sequence))
    for i, l in enumerate(open(options.xfilename).readlines()):
	entry = [float(x) for x in mysplit(l.strip(' \n,').replace(',',' '), ' ') ]
	for j, e in enumerate(entry):
	    bonuses[i][j] -=  e * adjust
    for i in range(len(bonuses)):
	bonusfile.write(','.join([str(x) for x in bonuses[i]]) + '\n')
    bonusfile.close()

    if options.shfilename:
	print 'Warning, both experimental bonus and shape constraint files specified.'
	print 'Going with the experimental bonus file and ignoring shape constraints.'
elif options.shfilename:
    bonusfilename ='/tmp/viennaRNA_bonuse_exp_for_shape.txt'
    bonusfile = open(bonusfilename, 'w')
    bonuses = []
    shbonuses = [0]*len(sequence)
    for i in range(len(sequence)):
	bonuses.append([0]*len(sequence))
    for l in open(options.shfilename).readlines():
	entry = [float(x) for x in mysplit(l.strip(), ' ')]
	entry[0] = int(entry[0])
	if entry[1] != -999:
	    for i in range(len(sequence)):
		bonuses[entry[0]-1][i] -= (options.sm*math.log(entry[1] + 1) + options.si) * adjust
		bonuses[i][entry[0]-1] -= (options.sm*math.log(entry[1] + 1) + options.si) * adjust
    for i in range(len(bonuses)):
	bonusfile.write(','.join([str(x) for x in bonuses[i]]) + '\n')
    bonusfile.close()
else:
    bonusfilename = ''
RNA.init_pf_fold(len(sequence))
if bonusfilename:
    print 'Adding bonuses'
    RNA.load_experimental_data(bonusfilename, len(sequence))
    RNA.load_pf_experimental_data(bonusfilename, len(sequence))
t = RNA.fold(sequence)
print 'MFE Structure: %s, dG: %s' % (t[0], t[1])
RNA.pf_fold(sequence)
m = []
for i in range(len(sequence)):
    m.append([0]*len(sequence))
for i in range(len(sequence)):
    for j in range(len(sequence)):
	if not math.isnan(RNA.get_pr(i+1,j+1)):
	    m[i][j] = RNA.get_pr(i+1,j+1)
	else:
	    m[i][j] = 0.0
for i in range(len(m)):
    outfile.write(','.join([str(x) for x in m[i]]) + '\n')
outfile.close()
