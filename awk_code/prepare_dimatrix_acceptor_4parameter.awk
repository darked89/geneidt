#!/usr/bin/gawk -f
#USAGE= $BIN/preparedimatrixacceptor4parameter.awk 18 19 20 $DATA/set2/sites/acceptor.geneid.dimatrix
BEGIN {
    preoffset = ARGV[1];
    newoffset = ARGV[2];
    posoffset = ARGV[3];
    
    ARGV[1] = ARGV[2] = ARGV[3] = "";
}
#  $2 !~ /A$/ #does not end in A
#  $2 !~ /^G/ #does not start with G
{
  if (($1 == preoffset && $2 !~ /A$/) || 
      ($1 == newoffset && $2 != "AG") || 
      ($1 == posoffset && $2 !~ /^G/))
      {
         print $1, $2, -9999;
      } 
  else 
	if (($1 == preoffset && $2   ~/A$/) ||
	    ($1 == newoffset && $2 == "AG"))
    print $1, $2, 0;
    else
       print;
}
