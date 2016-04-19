#!/usr/bin/env python2

"""

./cds2gff.awk ../pmar03.train.contig.gp.cds | less


# data structure
#$locus 1;
#$seqlen = 2;
#$strand = 3;
#$first_cds_exon_start = 4;

Purpose: takes a Golden Path annotation file and converts it to GFF
Problems(?): 
1) it is possible that the awk script is calculating phases of individual exons in a wrong way

2) the whole thing looks not sensible, since the coordinates should be for individual contigs. 
But then the names of the locus? are for corresponding to a single fasta...

3) no idea why but the frame here is always + ???

"""

import sys
from collections import deque

# gff defaults
source  =  "ANNOTATED"
GFF_SCORE = 1.000000
GFF_GROUP = 1


def calculate_exon_frame(exons_dict, current_exon_number, strand):
    pre_current_exon_len = 0
    total_num_exons = len(exons_dict)
    
    #print "###", current_exon_number
    if strand == "+":
		#assert current_exon_number > 1
		print "###", current_exon_number
		for exon_number in range(1, current_exon_number):
			tmp_start = exons_dict[exon_number]['start'] 
			tmp_end   = exons_dict[exon_number]['end'] 
			pre_current_exon_len = pre_current_exon_len + (tmp_end - tmp_start ) + 1
		exon_phase = pre_current_exon_len % 3
		
    elif strand == "-":
        #assert current_exon_number < total_num_exons
        for exon_number in range(current_exon_number+1, total_num_exons + 1):
            tmp_start = exons_dict[exon_number]['start'] 
            tmp_end   = exons_dict[exon_number]['end'] 
            pre_current_exon_len = pre_current_exon_len + (tmp_end - tmp_start ) + 1
        exon_phase = pre_current_exon_len % 3
    return exon_phase
    #print "###", current_exon_number, pre_current_exon_len, exon_phase
	  
cds_in_fn = sys.argv[1]

with open(cds_in_fn) as fh:
    counter = 1
    cds_in_dict = {}
    for line in fh:
		
		
		sl = line.split()
		seqname = sl[0]
		gene_label = "%s_%s" % (seqname, counter) 
		
		seqleng = sl[1]
		strand  = sl[2]
		exons_list = [int(x) for x in sl[3:]]
		exons_num  = len(exons_list)/2
		exons_dict = {}
		exons_queue = deque(exons_list)
		
		frame = 0
		reminder = 0
		
		for exon_counter in range(1, exons_num + 1):
			#list_counter = (exon_counter -1)
			exons_dict[exon_counter] = {}
			
			start = exons_queue.popleft()
			end   = exons_queue.popleft()
			
			exons_dict[exon_counter]['start'] = start #
			exons_dict[exon_counter]['end']   = end
			
			if exons_num == 1:
				exons_dict[exon_counter]['type']  = 'Single'
				exons_dict[exon_counter]['phase'] = 0
			elif exon_counter == 1:
				exons_dict[exon_counter]['type']   = 'First'
				exons_dict[exon_counter]['phase'] = 0 
				#remainder=($(fcds+1)-($fcds+frame)+1)%3;
				#nextframe=(3-($(fcds+1)-($fcds+frame)+1)%3)%3;
				  
				#remainder = (end - start + 1) % 3
				#nextframe = (3 - (end - (start + frame) + 1) % 3) %3
			elif exon_counter != exons_num:
				exons_dict[exon_counter]['type'] = 'Internal'
				exons_dict[exon_counter]['phase'] = calculate_exon_frame(exons_dict, exon_counter, strand)
				
			elif exon_counter == exons_num:
				exons_dict[exon_counter]['type'] = 'Last'
				exons_dict[exon_counter]['phase'] = calculate_exon_frame(exons_dict, exon_counter, strand)
			else:
				print "ERROR EXON", sl, exon_counter
			
			out_string = "%s\t%s\t%s\t%s\t%s" % (seqname, source, exons_dict[exon_counter]['type'], start, end) 
			out_string = "%s\t%s\t%s\t%s\t%s" % (out_string, GFF_SCORE, strand, exons_dict[exon_counter]['phase'], gene_label)
			#print "XXX", out_string
			print out_string
		#if exons_num != 1
				
		counter += 1
				
		#add some extra
		#print exons_num, exons_list
		#print exons_dict

		#f_cds   = sl[3]
''' 		
		if len(sl) == 6: #single exon gene
			printf gff_format, $locus, source, SINGLE, $fcds, $(fcds+1), SCORE, $strand, frame, $locus"_"group[$locus];  
			out_string = "%s\t%s\tSINGLE\t%s\t" % (seqname, source,
		
## deflt sourcrce name
#$source = 'ANNOTATED';



## gff formatted output
##$gff_format = "%s\t%s\t%s\t%d\t%d\t%f\t%s\t%d\t%s\n";


#lin while (<><>) {
    #chomp;	# strip record separator
    #@Fld = split(' ', $_, -1);

    #$fra = 0;
    #$groroup{$Fld[$locus]}++;

    ## check if single
   #(($#Fld+1) == $fcds + 1) {
	#printf $gffofof_format, $Fld[$locus], rce, $SINGLE, $Fld[$fcds], $],

	  #$Fld[$fcds + 1], $SCORE, $Fld[$strand $frame,

	  #$Fld[d[$locus] . '_' . $group{$Fld[$locus]};

	#next line
     #}

    #sumes fes first first exon, last last exon, othrs internal

    #int fir first exon
    #$reemainder = ($Fld[$fcds + 1] - ($Fld[$fcds] + $frame) + 1) % 3;
    #$nextframe = (3 - ($Fld[$fcds + 1] - ($Fld[$fcds] + $frame) + 1) % 3) % 3;

    #if d[$sd[$strand] eq '+') {
	#printf $gff_format, $Fld[$locus], $sou $FIRST, $Fld[$fcds],

	  #$Fld[$fcds + 1], $SCORSCORE, $Fld[$strand], $frame,

	  #$Fld[$locus] . '_' . $group{$Fld[$locus]};
    #}
    #else {
	#printf $gff_format, $Fld[$locus], $source, $TERMINAL, $Fld[$fcds],

	  #$Fld[$fcds + 1], $SCORE, $Fld[$strand], $remainder,

	  #$Fld[$locus] . '_' . $group{$Fld[$locus]};

	## prinnternernal exons
	#;
    #}
    #for ($i = $fcds + 2; $i < ($#Fld+1) - 2; $i += 2) 
	#$frame = $nextframe;
	#$remaemainder = ($Fld[$i + 1] - ($Fld[$i] + $frame) + 1) % 3;
	#$nextframe = (3 - ($Fld[$i + 1] - ($Fld[$i] + $frame) + 1) % 3) % 3;
	#if ($Fld[$strand] eq '+') {
	    #printf $gff_format, $Fld[$locus], $source, $INTERNAL, $Fld[$i],

	      #$Fld[$i + 1], $SCORE, $Fld[$strand], $frame,

	      #$Fld[$locus] . '_' . $group{$Fld[$locus]};
	#}
	#else {
	    #printf $gff_format, $Fld[$locus], $source, $INTERNAL, $Fld[$i],

	      #$Fld[$i + 1], $SCORE, $Fld[$strand], $remainder,

	      #$Fld[$locus] . '_' . $group{$Fld[$locus]};
	#}
    #}

   #rint ter terminal exon
    #$frame = $nextframe;
    #$remainder = ($Fld[$i + 1] - ($Fld[$i] + $frame) + 1) % 3;
    #if ($Fld[$strand] eq '+') {
	#printf $gff_format, $Fld[$locus], $source, $TERMINAL, $Fld[$i],

	  #$Fld[$i + 1], $SCORE, $Fld[$strand], $frame,

	  #$Fld[$locus] . '_' . $group{$Fld[$locus]};
    #}
    #else {
	#printf $gff_format, $Fld[$locus], $source, $FIRST, $Fld[$i],

	  #$Fld[$i + 1], $SCORE, $Fld[$strand], $remainder,

	  #$Fld[$locus] . '_' . $group{$Fld[$locus]};
    #}
#}
''' 

