#!/usr/bin/env python2

"""
replacement for frequency.awk

kmer-1/singe base only 2016.01.08
line to run orig:
gawk -f $path/frequency.awk 1 $true_seqs  > $true_seq_name.freq`;

"""
from __future__ import division

import sys

#kmer_len = sys.argv[1] is ever used something else?
kmer_len = 1 #hack at the moment


seq_dict   = {}
nucleotides_lst  = ('A', 'C', 'G', 'T', 'N')
nucleotides_clean  = ('A', 'C', 'G', 'T' ) 
#nucleotides_dict = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'ACGTsum': 0} 

kmer_num = 1 #hack at the moment
in_fn    = sys.argv[2]
#out_fn   = sys.argv[3]


old_seq_len = 0
counter     = 0
for line in open(in_fn).readlines():
    sl = line.split()
    sequence = sl[1]
    new_seq_len = len(sequence)
    if counter     != 0:
        pass
    else:
        old_seq_len = new_seq_len
    #elif new_seq_len != old_seq_len:
    #    print "seq len mismatch", line,
    #    break
    
    assert old_seq_len == new_seq_len
    assert new_seq_len > 0
    seq_dict[counter] = sequence
    counter += 1
    old_seq_len = new_seq_len

#print seq_dict

freq_dict = {}
#print new_seq_len
for base_pos in range(0, new_seq_len):
    freq_dict[base_pos] = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'ACGTsum': 0} 
    for key in seq_dict.keys():
         current_base = seq_dict[key][base_pos]
         assert current_base in nucleotides_lst
         if current_base != 'N':
            freq_dict[base_pos][current_base] += 1
            freq_dict[base_pos]['ACGTsum']    += 1

#nucleotides_lst_clean = nucleotides_lst - "N" 
#print "XXXX"

for position in freq_dict.keys():
    for base in nucleotides_clean:
        base_count       = freq_dict[position][base]
        total_base_count = freq_dict[position]['ACGTsum']
        base_frequency   = base_count / total_base_count
        #print "%s\t%s\t%s\t%.4f\t%s" % (base, position+1, freq_dict[position][base], freq_dict[position][base]/freq_dict[position]['ACGTsum'], freq_dict[position]['ACGTsum'])
        #print "%s\t%s\t%s\t%.4f" % (base, position, freq_dict[position][base], freq_dict[position][base]/freq_dict[position]['ACGTsum'])
        print "%s\t%s\t%s\t%.4f" % (position+1, base, base_frequency, base_count)

    #
    #print key, freq_dict[key]    
    
        




