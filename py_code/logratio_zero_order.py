#!/usr/bin/env python3

"""

#USAGE= gawk -f aux.awk Background_P_file, Observed_P_file
BEGIN{
    background_4cols_fn = ARGV[1]; 

    
    while (getline <background_4cols_fn > 0) { 
        bg_nucleotide = $1
        bg_positon    = $2
        bg_freq       = $4
        bg_freq_array[bg_nucleotide bg_position] = bg_freq
        }
      ## while(getline<ARGV[1]>0) # read background probabilities
      ##      BP[$1 $2]=$4;

    ARGV[1]="";
}
{   ##motif_4cols_fn = ARGV[2];
    
      if ($4!=0 && BP[$1 $2]!=0 && $4!=1) print $1, $2, log($4/BP[$1 $2]);
      if ($4!=0 && BP[$1 $2]!=0 && $4==1) print $1, $2, 0;
      if ($4==0 && BP[$1 $2]!=0) print $1, $2, -9999;
#else print $1, $2, 0;
}
"""
from __future__ import division
from math import log
import sys

nucleotides = 'ACGT'

background_4cols_fn = sys.argv[1]
signal_4cols_fn = sys.argv[2]

def read_4_cols(input_fn):
    tmp_frequency_dict = {}
    for line in open(input_fn).readlines():
        sl = line.split()
	### 202003 darked: 3 lines below differ from later  script
        position  = int(sl[0])
        base      = sl[1]
        frequency = float(sl[2])
        if position not in tmp_frequency_dict.keys():
            tmp_frequency_dict[position] = {}
        tmp_frequency_dict[position][base] = frequency
    return tmp_frequency_dict

#~ def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    #~ return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def isclose(a, b):
    tolerance=1e-09
    return abs(a-b) <= tolerance
    
background_freq_dict = read_4_cols(background_4cols_fn)
signal_freq_dict = read_4_cols(signal_4cols_fn)  

assert len(background_freq_dict) == len(signal_freq_dict)

for position in signal_freq_dict.keys():
    #~ ##DEBUG
    #~ if position in range(28,32):
        #~ print position, signal_freq_dict[position]
        #~ print position, background_freq_dict[position]
    for base in nucleotides:
        signal     = signal_freq_dict[position][base]
        background = background_freq_dict[position][base]

        assert background != 0.0

        if (isclose(signal, 0.0)):
            result = -9999
            
        elif (isclose(signal, 1.0)) and (not isclose(background, 1.0)):
            result = 0

        elif (not isclose(signal, 0.0)) and (not isclose(background, 1.0)):
            result = log(signal/background)
                        ## natural base e log
        else:
            print( "ERROR. bad values")
            print( position, base, signal_freq_dict[position][base])
            print( position, base, background_freq_dict[position][base])
            break
        #out_string = "%s\t%s\t%s" % (position, base, result)
        out_string = "%s %s %s" % (base, position,  result)
        print( out_string)

"""
      if ($4!=0 && BP[$1 $2]!=0 && $4!=1) print $1, $2, log($4/BP[$1 $2]);
      if ($4!=0 && BP[$1 $2]!=0 && $4==1) print $1, $2, 0;
      if ($4==0 && BP[$1 $2]!=0) print $1, $2, -9999;
            for base in nucleotides:
                tmp_frequency_dict[position] =
"""
