#!/usr/bin/env python2

"""
replacement for submatrix_order0.awk

#!/usr/bin/gawk -f

BEGIN {
    start = ARGV[1];
    end   = ARGV[2];
    
    ARGV[1] = ARGV[2] = "";
}
($2 >= start && $2 <= end ){
    $2=$2-start+1;
    print $2, $1, $3;
}

"""

import sys


start = int(sys.argv[1])
end   = int(sys.argv[2])
in_fn = sys.argv[3]

for line in open(in_fn).readlines():
	sl = line.split()
	col_two = int(sl[1])
	if (col_two >= start) and (col_two <= end):
		col_two = col_two - start + 1 
		print col_two, sl[0], sl[2]
