#!/usr/bin/gawk -f

# arguments: $temp_gp_cds_fn $temp_gp_gff 

BEGIN {
    while (getline < ARGV[1] >0) {
        len[$1]= $2;};
    
    OFS="\t";
    ARGV[1]= "";
    }

#main
{
    if (NR==1) {
        ant= $1;
        print  $1, $2, "Sequence", 1, len[ant], ".", ".", ".", "."
        };
    if ( $1!=ant) {
        #print  "# $ "; 
        ant= $1;
        print  $1, $2, "Sequence", 1, len[ant], ".", ".", ".", "."
        }; 
    print;
}

