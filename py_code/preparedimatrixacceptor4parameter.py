#!/usr/bin/env python2
'''

not sure if original awk produced best/intended matrix scores:
see also donor. Not consistent behaviour.
replicated original at this point

### org awk code below ###
BEGIN {
    st = ARGV[1];
    nd = ARGV[2];
    rd = ARGV[3];
    
    ARGV[1] = ARGV[2] = ARGV[3] = "";
}
#  $2 !~ /A$/ #does not end in A
#  $2 !~ /^G/ #does not start with G
{
  if (($1 == st && $2 !~ /A$/) || 
      ($1 == nd && $2 != "AG") || 
      ($1 == rd && $2 !~ /^G/))
      {
         print $1, $2, -9999;
      } 
  else 
	if (($1 == st && $2   ~/A$/) ||
	    ($1 == nd && $2 == "AG"))
    print $1, $2, 0;
    else
       print;
}
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
	if 	((col_one == pre_offset) and (seq[-1] != "A" )):
		print col_one, seq, -9999

	elif ((col_one == new_offset) and (seq != "AG")):
		print col_one, seq, -9999
	
	elif ((col_one == pos_offset) and (seq[0] != "G" )):
		print col_one, seq, -9999
    
	elif ((col_one == pre_offset) and  (seq[-1] == "A")):
		print col_one, seq, 0
	elif ((col_one == new_offset) and  (seq     == "AG")):
		print col_one, seq, 0
	#
	#elif ((col_one == pos_offset ) and (seq[0] == "G")):
	#	print col_one, seq, 0
	else:
		print col_one, seq, col_three


