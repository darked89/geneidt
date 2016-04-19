#!/usr/bin/env python2
"""
replacement of corresponding GAWK

differences:
1) tabs instead of spaces?
2) no extra "," after the bunch of exons
3) no margins for true exons (first/last) 

"""

import sys

MARGIN = 400 # context can be changed according to how much sequence we want on 3' and 5' of gene
input_fn = sys.argv[1]

genes_dict = {}

for line in open(input_fn).readlines():
	sl = line.split()
	gene_name = sl[8]
	if gene_name not in genes_dict.keys():
		genes_dict[gene_name] = {}
		genes_dict[gene_name]['num_exons'] = 0
		genes_dict[gene_name]['exons']     = {}
		genes_dict[gene_name]['exons'][1]  = [] #list of a start and end
	
	#missing checks for gff file sanity
	#like same gene name with exons on different contigs
	#sorted gff or not? 
	genes_dict[gene_name]['contig_id']  = sl[0]
	genes_dict[gene_name]['strand']     = sl[6]
	genes_dict[gene_name]['num_exons'] += 1
	current_exon_number = genes_dict[gene_name]['num_exons']
	genes_dict[gene_name]['exons'][current_exon_number] = [int(sl[3]), int(sl[4])] #exon start plus exon end


for gene_name in genes_dict.keys():
	num_exons = genes_dict[gene_name]['num_exons']
	gp_string = "%s\t%s\t%s\t" % (gene_name, genes_dict[gene_name]['contig_id'],genes_dict[gene_name]['strand'] )
	 
	gp_positions = "%s\t%s\t%s\t%s\t%s\t" % (\
	genes_dict[gene_name]['exons'][1][0] - (MARGIN+1), \
	genes_dict[gene_name]['exons'][num_exons][1] + (MARGIN+1), \
	genes_dict[gene_name]['exons'][1][0] - 1, \
	genes_dict[gene_name]['exons'][num_exons][1], \
	num_exons)
	#gp
	
	exons_start_positions = ""
	exons_end_positions   = ""
	for exon_number  in range(1, num_exons+1):
		if exon_number == 1:
			exons_start_positions = genes_dict[gene_name]['exons'][exon_number][0] #- (MARGIN+1)
			exons_end_positions   = genes_dict[gene_name]['exons'][exon_number][1] 
		else:
			exons_start_positions = "%s,%s"  %  (exons_start_positions, genes_dict[gene_name]['exons'][exon_number][0])
			exons_end_positions   = "%s,%s"  %  (exons_end_positions,   genes_dict[gene_name]['exons'][exon_number][1])
	
	#print gene_name, genes_dict[gene_name]
	print gp_string, gp_positions, exons_start_positions, exons_end_positions
	
	
	
