#!/usr/bin/gawk -f 

BEGIN {
	while (getline < ARGV[1] > 0) 
		{
		length[$1] = $2;
		};

	ARGV[1] = "";
	OFS     = "\t";
}


{
if (NR == 1) 
{
	ant = $1;
	print $1, $2, "Sequence", 1, length[ant], ".", ".", ".", "." 
};

if ( $1 != ant ) 
	{
	print "\#\$";
	ant = $1;
	print $1, $2, "Sequence", 1, length[ant], ".", ".", ".", "."
	};
	print 
}' 

#$tempcdsgp $tempgp_gff |";

