#!/usr/bin/gawk -f

# invoked as:
# 2562-        my $gpcontig = "";
# 2563-        open LOCID,
# 2564:"$path/multiple_annot2one.awk species=$species leng=$lengp $tempcdsgp |"
# 2565-          ;    #possible error dk 2014.11.08
# 2566-        while (<LOCID>) { $gpcontig .= $_; }
# 2567-        close LOCID;

#
# converts geneid gff output into DAS mysql compliance data

#dk: pointless, single space is the default out fiels separator?
#BEGIN{
#	OFS=" ";
#    }

{


#NR stands for number of records in file (lines?)
if (NR == 1) {
	tot = split($0, array, " ");
	seq = $1;
	add = $2;
	for (coord=4; coord <= tot; coord++)
		{
		 s = s""array[coord]" "
		 };

tot = 0;
print species, leng, $3, s;
s = "";
oldadd = 0;
}


if ($1 != seq) {
    tot = split($0, array," ");
	seq = $1;
	for (coord=4; coord<=tot; coord++) {
		val = array[coord] + add;
		s = s""val" "
		};
	
	print species, leng, $3, s;
	s = "";
	add = (add + $2);
	tot = 0;
	}
 
}
