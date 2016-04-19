#!/usr/bin/env python2
'''

#USAGE:

./preparedimatrixdonor4parameter.py 3 4 5 canonical.donor.matrix.pre

# $pre_offset $new_offset $pos_offset

temp solution. 
not sure if original awk produced best/intended matrix scores:
3rd di-nucleotide terminating in G(first G in GT donor) score 0?
5th di-nucleotide starting with  T (last T in GT donor) score 0?
replicated original at this point

'''
import sys

pre_offset = int(sys.argv[1])
new_offset = int(sys.argv[2])
pos_offset = int(sys.argv[3])
in_fn     = sys.argv[4]

#print pre_offset, new_offset, pos_offset, in_fn

for line in open(in_fn).readlines():
	sl = line.split()
	col_one   = int(sl[0])
	seq       = sl[1]
	col_three = float(sl[2])
	 
	#print sl
	#print pre_offset, new_offset, pos_offset, 
	#print col_one, seq, seq[-1], col_three
	if 	((col_one == pre_offset) and (seq[-1] != "G" )):
		print col_one, seq, -9999

	elif ((col_one == new_offset) and (seq != "GT")):
		print col_one, seq, -9999
	
	elif ((col_one == pos_offset) and (seq[0] != "T" )):
		print col_one, seq, -9999
    
	#elif ((col_one == pre_offset) and  (seq[-1] == "G")):
	#	print col_one, seq, 0
	elif ((col_one == new_offset) and  (seq     == "GT")):
		print col_one, seq, 0
	#
	#elif ((col_one == pos_offset ) and (seq[0] == "T")):
	#	print col_one, seq, 0
	else:
		print col_one, seq, col_three


