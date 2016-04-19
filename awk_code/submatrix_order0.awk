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
