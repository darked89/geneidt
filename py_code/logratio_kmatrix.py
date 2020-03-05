#!/usr/bin/env python3

"""
#USAGE= gawk -f aux.awk Background_P_file, Observed_P_file
BEGIN{
  while(getline<ARGV[1]>0) # read background probabilities
    BP[$1 $2]=$3;

  ARGV[1]="";
}
{
  if ($3!=0 && BP[$1 $2]!=0 ) print $1, $2, log($3/BP[$1 $2]);
  else print $1, $2, 0;
}
"""

from __future__ import division
from math import log
import sys
from itertools import product

nucleotides = 'ACGT'

background_4cols_fn = sys.argv[1]
signal_4cols_fn = sys.argv[2]

def read_3_cols(input_fn):
    tmp_frequency_dict = {}
    result_freq_dict = {}
    for line in open(input_fn).readlines():
        sl = line.split()
        position  = int(sl[0])
        kmer      = sl[1]
        frequency = float(sl[2])
        if position not in tmp_frequency_dict.keys():
            tmp_frequency_dict[position] = {}
        tmp_frequency_dict[position][kmer] = frequency
    ##result_freq_dict 
    return tmp_frequency_dict

#~ def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    #~ return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def isclose(a, b):
    tolerance=1e-09
    return abs(a-b) <= tolerance

four_bases = 'ACGT'
kmer_plus_list2 = []
for element in product(four_bases, repeat=2):
    my_string = ''.join(x for x in element)
    #print my_string
    kmer_plus_list2.append(my_string)
    
background_freq_dict = read_3_cols(background_4cols_fn)
signal_freq_dict = read_3_cols(signal_4cols_fn)  

assert len(background_freq_dict) == len(signal_freq_dict)

for position in signal_freq_dict.keys():
    #~ ##DEBUG
    #~ if position in range(28,32):
        #~ print position, signal_freq_dict[position]
        #~ print position, background_freq_dict[position]
    for kmer in kmer_plus_list2:
        signal     = signal_freq_dict[position][kmer]
        background = background_freq_dict[position][kmer]

        assert background != 0.0

        if isclose(signal, 0.0) or isclose(background, 0.0):
            result = 0
        

        elif (not isclose(signal, 0.0)) and (not isclose(background, 1.0)):
            result = log(signal/background)
                        ## natural base e log
        else:
            print( "ERROR. bad values")
            print( position, kmer, signal_freq_dict[position][kmer])
            print( position, kmer, background_freq_dict[position][kmer])
            break
        out_string = "%s\t%s\t%s" % (position, kmer, result)
        out_string = "%s %s %s" % (position, kmer, result)
        #out_string = "%s %s %s" % (kmer, position,  result)
        print( out_string )
