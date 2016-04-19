#!/usr/bin/env python2
'''

not sure if original awk produced best/intended matrix scores:
see also donor. Not consistent behaviour.
replicated original at this point

### org awk code below ###

#USAGE= $BIN/preparetrimatrixstart4parameter.awk 4 5 6 7 $DATA/set2/sites/starts.geneid.dimatrix
BEGIN {
    st = ARGV[1];
    nd = ARGV[2];
    rd = ARGV[3];
    th = ARGV[4];   
    ARGV[1] = ARGV[2] = ARGV[3] = ARGV[4] = "";
}
{
  if (($1 == st && $2 !~ /A$/) || 
      ($1 == nd && $2 != "AT") || 
      ($1 == rd && $2 != "TG") || 
      ($1 == th && $2 !~ /^G/))
      {
         print $1, $2, -9999;
    } else 
    if (($1 == st && $2 ~  /A$/) ||\
        ($1 == nd && $2 == "AT") || 
        ($1 == rd && $2 == "TG"))
    print $1, $2, 0;
  else
print;
}
'''

import sys

A_offset = int(sys.argv[1])
B_offset = int(sys.argv[2])
C_offset = int(sys.argv[3])
D_offset = int(sys.argv[4])

in_fn     = sys.argv[5]

#print pre_offset, new_offset, pos_offset, in_fn

for line in open(in_fn).readlines():
	sl = line.split()
	col_one   = int(sl[0])
	seq       = sl[1]
	col_three = float(sl[2])
	 
	#print sl
	#print pre_offset, new_offset, pos_offset, 
	#print col_one, seq, seq[-1], col_three
	#-9999 vals
	if 	((col_one == A_offset) and (seq[-1] != "A" )):
		print col_one, seq, -9999

	elif ((col_one == B_offset) and (seq != "AT")):
		print col_one, seq, -9999
	
	elif ((col_one == C_offset) and (seq != "TG" )):
		print col_one, seq, -9999

    elif ((col_one == D_offset) and (seq[0] != "G" )):
		print col_one, seq, -9999
	
	#zeros	
	elif ((col_one == A_offset) and  (seq[-1] == "A")):
		print col_one, seq, 0
	elif ((col_one == B_offset) and  (seq     == "AT")):
		print col_one, seq, 0
	
	elif ((col_one == C_offset) and  (seq     == "TG")):
		print col_one, seq, 0

	#
	#elif ((col_one == D_offset ) and (seq[0] == "G")):
	#	print col_one, seq, 0
	else:
		print col_one, seq, col_three


