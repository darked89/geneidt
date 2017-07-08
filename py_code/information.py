#!/usr/bin/env python2

from __future__ import division

import sys

from math import log


true_freq_fn = sys.argv[1]
bckg_freq_fn = sys.argv[2]



#print "## from information.py: "
#print "## input args1: bckg_freq_fn", bckg_freq_fn
#print "## input arg12: true_freq_fn", true_freq_fn 


trueseq_dict = {}
bckgseq_dict = {}

counter    = 0
kmer_size  = 0
for line in open(true_freq_fn).readlines():
    sl = line.split()
    position = int(sl[0])
    kmer     = sl[1]
    if position not in  trueseq_dict.keys():
        trueseq_dict[position] = {}
    trueseq_dict[position][kmer] = float(sl[2])
    #my_key = "%s_%s" % (sl[0],sl[1])
    #trueseq_dict[my_key] = float(sl[3])
    if counter == 0:
        kmer_size = len(sl[0])
    counter += 1

#debug line, not needed in the output
#print counter
#print trueseq_dict


new_counter = 0
for line in open(bckg_freq_fn).readlines():
    if new_counter < counter:
        sl = line.split()
        position = int(sl[0])
        kmer     = sl[1]
        if position not in bckgseq_dict.keys():
            bckgseq_dict[position] = {}
        bckgseq_dict[position][kmer] = float(sl[2])
    
        #my_key = "%s_%s" % (sl[0],sl[1])
        #bckgseq_dict[my_key] = float(sl[3])
        new_counter += 1

#print bckgseq_dict
#debug line, not needed in the output
#print new_counter, len(bckgseq_dict), len(trueseq_dict), kmer_size

#make sure that somehow we are not iterating beyond the background kmer length
assert len(trueseq_dict.keys()) <=  len(bckgseq_dict.keys())

for position in trueseq_dict.keys():
    c_variable = 0.0;
    for kmer in trueseq_dict[position].keys():
        #print "XXX", kmer
        true_freq = trueseq_dict[position][kmer]
        bckg_freq = bckgseq_dict[position][kmer]
        #if bckg_freq <= 0:
        #    print position, kmer, bckg_freq[position][kmer]
        try:
            temp_info_val = true_freq*log(true_freq/bckg_freq)/log(2)
            #Sprint "Y", position, kmer, temp_info_val
        except: 
            temp_info_val = 0 #not sure
        c_variable += temp_info_val
        #print kmer, position,  bckg_freq, true_freq, temp_info_val
        #debug line, not needed in the output
        #print kmer, position, temp_info_val
    print position, c_variable
    c_variable = 0.0
'''    
        true_freq = trueseq_dict[kmer_position]
    if bckg_freq != 0: #not sure....
        temp_info = true_freq*log(true_freq/bckg_freq)/log(2)
    else:
        temp_info = 0
    
    if ((bckg_freq == 0) or (true_freq == 0)):
        c_variable = 
'''            
