#!/usr/bin/gawk -f
#USAGE= $BIN/preparedimatrixdonor4parameter.awk 3 4 5 $DATA/set2/sites/donor.geneid.dimatrix
# $preoffset $newoffset $posoffset
BEGIN {
    preoffset = ARGV[1];
    newoffset = ARGV[2];
    posoffset = ARGV[3];
    ARGV[1] = ARGV[2] = ARGV[3] = "";
}
{
  if (($1 == preoffset  && $2 !~ /G$/) ||
      ($1 == newoffset  && $2 != "GT") || 
      ($1 == posoffset  && $2 !~ /^T/) ) 
      {
         print $1, $2, -9999;
    } 
   else 
    if (($1 == preoffset && $2 ~ / G$/) || 
        ($1 == newoffset && $2 == "GT") )
    print $1, $2, 0;
  else
print;
}
