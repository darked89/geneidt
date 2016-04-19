#!/usr/bin/env python2

"""
replacement for submatrix.awk

"""

import sys


start = int(sys.argv[1])
end   = int(sys.argv[2])
in_fn = sys.argv[3]

for line in open(in_fn).readlines():
	sl = line.split()
	col_one = int(sl[0])
	if (col_one >= start) and (col_one <= end):
		col_one = col_one - start + 1 
		print col_one, sl[1], sl[2]
