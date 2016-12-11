#!/usr/bin/env python2

"""
simple replacement for gi_trainer for geneid

ideas to implement:
1) gff parse taken partialy from
https://gist.github.com/slowkow/8101481

2) interval/gene sanity checking:
https://www.biostars.org/p/99/#101
https://bitbucket.org/james_taylor/bx-python/wiki/Home


"""

import sys, os
import re
#import yaml

from pyfaidx import Fasta



#~ conf_yaml_fn = sys.argv[1]

#~ with open(conf_yaml_fn, 'r') as ymlfile:
    #~ cfg = yaml.load(ymlfile)


#~ #print cfg


def check_fasta(fasta_fn):
    genome_fas = Fasta(fasta_fn , sequence_always_upper=True)
    genome_fas_dict = {}
    
    for contig_name in genome_fas.keys():
        #print contig_name, len(genome_fas[contig_name][:])
        genome_fas_dict[contig_name] = len(genome_fas[contig_name][:])
        #print "%key"
    return len(genome_fas.keys())


#~ #temp_gff_list     = []

#~ gff_values_counts = {}
#~ col3_vals = cfg['gff']['col3_vals']

#~ for value in col3_vals:
    #~ gff_values_counts[value] = 0    

#~ def gff_reader(gff_fn):
    
    #~ #temp_gff_list     = []
    #~ transcript_set = set()
    #~ transcript_dict = {}
    #~ gff_values_counts = {}
    #~ int_values =  cfg['gff']['int_values']
    
    #~ for value in col3_vals:
        #~ gff_values_counts[value] = 0
    #~ temp_transcript_id = ''    
    
    #~ for i, gffline in enumerate(gfflines(gff_fn)):
        #~ for int_value in int_values:
            #~ gffline[int_value] = int(gffline[int_value])
            
        #~ gffline['transcript_id']
        #~ gff_values_counts[gffline['feature']] += 1
        
        #~ temp_transcript_id = gffline['transcript_id']
        
        #~ if temp_transcript_id not in transcript_dict.keys():
            #~ transcript_dict[temp_transcript_id] = []
        
        #~ transcript_dict[temp_transcript_id].append(gffline)
        
        #~ transcript_set.add(gffline['transcript_id'])
        
        #~ #temp_gff_list.append(gffline)    
    
    #~ print "processed %s lines in %s " % (i, gff_fn)
    #~ print "obtained %s transcript names" % (len(transcript_set))
    #~ print "values/exon types counts: "
    #~ print gff_values_counts
    
    #~ exon_numbers       = {}
    #~ plus_strand_trans  = 0
    #~ minus_strand_trans = 0
    #~ for transcript_id in transcript_dict.keys():
        #~ current_record = transcript_dict[transcript_id]
        #~ num_of_exons = len(current_record)
        #~ #print num_of_exons
        #~ if num_of_exons not in exon_numbers.keys():
            #~ exon_numbers[num_of_exons] = 1
        #~ else:
            #~ exon_numbers[num_of_exons] += 1
        
        #~ if  current_record[0]['strand'] == '+':
            #~ plus_strand_trans += 1
        
        #~ elif current_record[0]['strand'] == '-':
            #~ minus_strand_trans += 1
        #~ else:
            #~ print "bad record, no strand: ", current_record
    
    #~ print "plus_strand_trans: %s " % (plus_strand_trans)
    #~ print "minus_strand_trans: %s " % (minus_strand_trans)
    
    #~ for key in exon_numbers.keys():
        #~ print "exon transcriots",  key, exon_numbers[key] 
    
    
    #~ return transcript_dict
    
    
    
    



#~ def gfflines(filename):
    #~ """Open an optionally gzipped GTF file and generate a dict for each line.
    #~ """
    #~ #fn_open = gzip.open if filename.endswith('.gz') else open

    #~ with open(filename) as fh:
        #~ for line in fh:
            #~ if line.startswith('#'):
                #~ continue
            #~ else:
                #~ yield parse(line)
            
                

#~ def parse(line):
    #~ """Parse a single GTF line and return a dict.
    #~ """
    #~ GTF_HEADER = cfg['GTF_HEADER']
    #~ result = {}

    #~ fields = line.rstrip().split('\t')
    #~ if len(fields) !=9:
        #~ print line,
        #~ print fields
    #~ for i, col in enumerate(GTF_HEADER):
        #~ result[col] = _get_value(fields[i])
    #~ return result
    
#~ def _get_value(value):
    #~ if not value:
        #~ return None

    #~ # Strip double and single quotes.
    #~ value = value.strip('"\'')

    #~ # Return a list if the value has a comma.
    #~ if ',' in value:
        #~ value = re.split(R_COMMA, value)
    #~ # These values are equivalent to None.
    #~ elif value in ['', '.', 'NA']:
        #~ return None

    #~ return value
    
    

#~ gff_fn = cfg['gff']['genes_fn']
#~ transcript_dict = gff_reader(gff_fn)

#~ print "transcript_dict size:", len(transcript_dict)


#~ genome_fas = Fasta( cfg['genome_fasta'], sequence_always_upper=True)
#~ genome_fas_dict = {}
#~ for key in genome_fas.keys():
    #~ genome_fas_dict[key] = len(genome_fas[key][:])

#~ #print "got here"
#~ #print  '10274.m000709', genome_fas_dict['10274.m000709']
#~ #extract sequences

#~ def extract_splices(transcript_dict, genome_fas, genome_fas_dict, splice_type):

    #~ def plus_strand_feat(chrom, position, offset):
        #~ chrom_size = genome_fas_dict[chrom]
        #~ #print chrom, chrom_size, position, offset
        #~ prefix_len  = 0
        #~ postfix_len = 0
        #~ prefix_str  = ""
        #~ postfix_str = ""
        
        #~ if (position  - offset)  < 1:
            #~ prefix_len = abs(position - offset) + 1
            #~ prefix_str = "N" * prefix_len
            #~ temp_seq = prefix_str + genome_fas[chrom][:position + offset +1]
        #~ elif position+offset > chrom_size:
            #~ postfix_len = (position+offset) - chrom_size
            #~ postfix_str = "N" * postfix_len
            #~ temp_seq = genome_fas[chrom][position - offset:chrom_size] + postfix_str
        #~ else:
            #~ temp_seq = genome_fas[chrom][position - offset-1:position + offset-1]
        #~ return temp_seq, len(temp_seq)
        
    #~ def minus_strand_feat(chrom, position, offset):
        #~ chrom_size = genome_fas_dict[chrom]
        #~ #print chrom, chrom_size, position, offset
        #~ prefix_len  = 0
        #~ postfix_len = 0
        #~ prefix_str  = ""
        #~ postfix_str = ""
        
        #~ if (position  - offset)  < 1:
            #~ prefix_len = abs(position-offset) + 1 
            #~ prefix_str = "N" * prefix_len
            #~ temp_seq = prefix_str + genome_fas[chrom][:position+offset].complement
        #~ elif position+offset > chrom_size:
            #~ postfix_len = (position+offset) - chrom_size
            #~ postfix_str = "N" * postfix_len
            #~ temp_seq = genome_fas[chrom][position-offset:chrom_size].complement + postfix_str
        #~ else:
            #~ temp_seq = genome_fas[chrom][position-offset:position+offset].complement
        #~ temp_seq = temp_seq[::-1] #reversing
        #~ return temp_seq, len(temp_seq)
    
    
    #~ #counter = 0
    #~ def get_acceptors(transcript_list):
        #~ if transcript_list == None:
            #~ transcript_list = transcript_dict.keys()
        #~ else:
            #~ print transcript_list
     
        #~ for transcript_id in transcript_list:
            #~ current_record = transcript_dict[transcript_id]
            #~ num_of_exons = len(current_record)
            #~ #print transcript_id, num_of_exons
            
            #~ if num_of_exons > 1:
                #~ strand = current_record[0]['strand']
                #~ if strand == '+':
                    #~ for exon_num in range(1, num_of_exons):
                        #~ chrom    = current_record[exon_num]['seqname']
                        #~ position = current_record[exon_num]['start']
                        #~ seq, seq_len = plus_strand_feat(chrom, position, offset)
                        #~ print "%s_Fex%s_%s_%s_%s\t%s" % (transcript_id, exon_num, label, seq_len, seq[offset-2:offset], seq)
                #~ elif strand == '-':
                    #~ for exon_num in range(0, num_of_exons-1):
                        #~ chrom    = current_record[exon_num]['seqname']
                        #~ position = current_record[exon_num]['end']
                        #~ seq, seq_len = minus_strand_feat(chrom, position, offset)
                        #~ print "%s_Rex%s_%s_%s_%s\t%s" % (transcript_id, exon_num, label, seq_len, seq[offset-2:offset], seq)
                #~ else:
                    #~ print "XXX strand error: no strand info for transcript_id"
                    #~ print "XXX", current_record
    #~ def get_donors(transcript_list):
        #~ if transcript_list == None:
            #~ transcript_list = transcript_dict.keys()
     
        #~ for transcript_id in transcript_list:
            #~ current_record = transcript_dict[transcript_id]
            #~ num_of_exons = len(current_record)
            #~ #print transcript_id, num_of_exons
            
            #~ if num_of_exons > 1:
                #~ strand = current_record[0]['strand']
                #~ if strand == '+':
                    #~ for exon_num in range(0, num_of_exons-1):
                        #~ chrom    = current_record[exon_num]['seqname']
                        #~ position = current_record[exon_num]['end']
                        #~ seq, seq_len = plus_strand_feat(chrom, position, offset)
                        #~ print "%s_Fex%s_%s_%s_%s\t%s" % (transcript_id, exon_num, label, seq_len, seq[offset+1:offset+3], seq)
                #~ elif strand == '-':
                    #~ for exon_num in range(1, num_of_exons-1):
                        #~ chrom    = current_record[exon_num]['seqname']
                        #~ position = current_record[exon_num]['start']
                        #~ seq, seq_len = minus_strand_feat(chrom, position, offset)
                        #~ print "%s_Rex%s_%s_%s_%s\t%s" % (transcript_id, exon_num, label, seq_len, seq[offset+1:offset+3], seq)
                #~ else:
                    #~ print "XXX strand error: no strand info for transcript_id"
                    #~ print "XXX", current_record            
        #~ #the AG from
    
    #~ #executing the functions
    #~ offset   = 30
    #~ position = 0
    #~ label = splice_type[:3]
    #~ transcript_list = transcript_dict.keys()
    
    #~ if  splice_type == "acceptors":
        #~ #AG types, exon numbering  from 0
        #~ #assuming position sorted not nested genes
                #~ #F_strand_exons_num_range = ('start', 1, "last")
        #~ #R_strand_exons_num_range = ('end',   0, "second_last")
        #~ get_acceptors(transcript_list)
        
    #~ if  splice_type == "donors":
        #~ #GT type
        #~ #F_strand_exons_num_range = ('end',   0, "second_last")
        #~ #R_strand_exons_num_range = ('start', 1, "last")
        #~ get_donors(transcript_list)
        

#~ def extract_cds_seqs(*transcript_list):
    #~ transcript_list = transcript_dict.keys()
    #~ #debug
    #~ #transcript_list = ['model.1879.m002022', 'foobar']
    #~ #pattern for longest ORF in transcript
    #~ orf_finder = re.compile(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)')
    
    #~ def get_exon_seq(exon_record, strand):            
        #~ current_exon = exon_record
        #~ chrom    = current_exon['seqname']
        #~ start    = current_exon['start']
        #~ end      = current_exon['end']
        #~ if strand == "+":                
            #~ exon_seq = genome_fas[chrom][start-1:end]
        #~ if strand == "-":
            #~ exon_seq = genome_fas[chrom][start-1:end].complement
            #~ exon_seq = exon_seq[::-1]
        #~ #print "RRR", len(exon_seq)
        #~ #print "RRR\t%s" % (exon_seq)
        #~ return exon_seq
        
    
    
    #~ for transcript_id in transcript_list:
        #~ current_record = transcript_dict[transcript_id]
        #~ num_of_exons = len(current_record)
        #~ strand = current_record[0]['strand']
        
        #~ transcript_seq = ""
        #~ cds_seq        = "" 
        #~ if num_of_exons == 1:
            #~ current_exon = current_record[0]
            #~ transcript_seq = "%s" % get_exon_seq(current_exon, strand)
        #~ else:
            #~ for exon in current_record:
                #~ #print "$$$$", exon
                #~ transcript_seq += "%s" % get_exon_seq(exon, strand)
        #~ print ">%s_from_%s_exons_%s" % (transcript_id, num_of_exons, strand)
        #~ print "%s" % transcript_seq
        #~ try:        
            #~ cds_seq = max(orf_finder.findall(transcript_seq), key = len)
            #~ print "cds_seq", cds_seq
        #~ except:
            #~ print "error: transcript_id", transcript_id

#~ extract_splices(transcript_dict, genome_fas, genome_fas_dict, "acceptors")
#~ extract_splices(transcript_dict, genome_fas, genome_fas_dict, "donors")
#~ extract_cds_seqs()
