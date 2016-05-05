#!/usr/bin/env python2

import sys, os

from random import randint

from pyfaidx import Fasta

kmer_size = 60

fasta_in = sys.argv[1]

genome_fas = Fasta( fasta_in, sequence_always_upper=True)
genome_fas_dict = {}
cumulative_size_list = []
contigs         = []
cumulative_lenght = 0 
#total_bp_genomic = 0

for key in genome_fas.keys():
    contig_length = len(genome_fas[key][:])
    cumulative_lenght += contig_length
    cumulative_size_list.append(cumulative_lenght)
    contigs.append(key)
    genome_fas_dict[key] = len(genome_fas[key][:])


total_bp_genomic = cumulative_lenght
#print cumulative_lenght
#print cumulative_size_list
#print contigs

#previous_contig = ""
#previous_contig_end_pos = 0
#current_contig_end_pos  = 0 
#contig_index   = 0
kmer_num_counter = 0
while kmer_num_counter < 500000:
#while kmer_num_counter < 5:
    #current_end_pos  = 0
    #previous_contig = ""
    rand_start_pos = randint(1,total_bp_genomic - kmer_size)
    #print rand_start_pos
    for current_end_pos in cumulative_size_list:
        if rand_start_pos  > (current_end_pos - kmer_size):
            pass
        else:

            contig_idx = cumulative_size_list.index(current_end_pos)
            #print contig_idx
            contig_name = contigs[contig_idx]
            #print contig_name
            if contig_idx != 0:
                previous_cumulative_end = cumulative_size_list[contig_idx - 1]
            else:
                previous_cumulative_end = 0
            start_in_contig = rand_start_pos - previous_cumulative_end
            end_in_contig   =  start_in_contig +  kmer_size
            #print rand_start_pos, previous_cumulative_end, start_in_contig, end_in_contig
            try:
                kmer_seq = genome_fas[contig_name][start_in_contig:end_in_contig]
                tmp_out_str = "kmer%s\t%s\t%s" %(kmer_num_counter, kmer_seq, len(kmer_seq))
                tmp_out_str = "km%s\t%s" %(kmer_num_counter, kmer_seq)
                print tmp_out_str
                #print kmer_seq
            except:
                kmer_num_counter -= 1
                
            #
            #print tmp_out_str
            kmer_num_counter += 1
            break


'''           
    while rand_start_pos  > current_end_pos - kmer_size:
        current_contig = contigs[contig_index]
        cur_end_pos  = cumulative_size[contig_index] 
        contig_index += 1
    selected_contig = current_contig
    start_in_selected_contig = rand_start_pos - cumulative_size[contig_index-1]
    
    tmp_out_str = "5s\t%s" %(kmer_num_counter, kmer_seq)
    print tmp_out_str
''' 

    

    

