#!/usr/bin/python

import RNA

##############################################
# Run RNA fold
# and sum over base pair probability matrix
# to get probability of each base i forming
# a base pair with anything else.
##############################################
def get_base_exposure( sequence, temperature ):

    RNA.cvar.temperature = temperature
    RNA.init_pf_fold( 1 )
    RNA.pf_fold( sequence )
    base_exposure = []

    for i in range( len( sequence ) ):

        base_pair_prob = 0.0

        for j in range( i-1 ):
            prob = RNA.get_pr( i+1, j+1 )
            base_pair_prob += prob
        for j in range( i, len(sequence) ):
            prob = RNA.get_pr( i+1, j+1 )
            base_pair_prob += prob

        base_exposure.append( 1.0 - base_pair_prob )

    return base_exposure


def get_tapestry( sequence, mutpos = [], offset = 0, temperature = 24.0, only_library1 = 0 ):

    tapestry  = []
    mutations = []

    if len( mutpos ) == 0:
        mutpos = range( len( sequence ) )

    exposure = get_base_exposure( sequence, temperature )
    tapestry.append( exposure )
    mutations.append( 'WT' )

    sequence_new = ''
    for i in range( len( sequence) ):
        if sequence[i] == 'T':
            sequence_new += 'U'
        else:
            sequence_new += sequence[i]
    sequence = sequence_new

    # 3 possible mutations
    mutation_keys = []
    mutation_keys.append( {'A':'U', 'U':'A', 'G':'C', 'C':'G' } )
    mutation_keys.append( {'A':'C', 'U':'C', 'G':'A', 'C':'A' } )
    mutation_keys.append( {'A':'G', 'U':'G', 'G':'U', 'C':'U' } )

    num_libraries = 3
    if only_library1: num_libraries = 1

    for library in range( num_libraries ):

        for i in range( len( sequence ) ):

            if ( (i + 1 + offset) not in mutpos ): continue

            sequence_mut = sequence[:i] + mutation_keys[library][ sequence[i] ] + sequence[(i+1):]
            #print sequence_mut
            exposure = get_base_exposure( sequence_mut, temperature )


            tapestry.append( exposure )
            mutation = '%s%d%s' % (sequence[i], i+1+offset, mutation_keys[library][sequence[i]] )
            #print mutation
            mutations.append( mutation )

            #for m in range( len( exposure ) ):
            #    if exposure[m ] < 0.0: print "PROBLEM!! ", mutation, exposure[ m ]

    return (tapestry, mutations)
