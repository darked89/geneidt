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
import yaml

from pyfaidx import Fasta

from Bio.Seq import Seq

#~ from string import maketrans
#~ fw = 'ACGTN'
#~ rv = 'TGCAN'
#~ nucl_table = maketrans(fw, rv)

conf_yaml_fn = sys.argv[1]

with open(conf_yaml_fn, 'r') as ymlfile:
    cfg = yaml.load(ymlfile)


#print cfg



#temp_gff_list     = []

gff_values_counts = {}
col3_vals = cfg['gff']['col3_vals']

for value in col3_vals:
    gff_values_counts[value] = 0    

def gff_reader(gff_fn):
    
    #temp_gff_list     = []
    transcript_set = set()
    transcript_dict = {}
    gff_values_counts = {}
    int_values =  cfg['gff']['int_values']
    
    for value in col3_vals:
        gff_values_counts[value] = 0
    temp_transcript_id = ''    
    
    for i, gffline in enumerate(gfflines(gff_fn)):
        for int_value in int_values:
            gffline[int_value] = int(gffline[int_value])
            
        gffline['transcript_id']
        gff_values_counts[gffline['feature']] += 1
        
        temp_transcript_id = gffline['transcript_id']
        
        if temp_transcript_id not in transcript_dict.keys():
            transcript_dict[temp_transcript_id] = []
        
        transcript_dict[temp_transcript_id].append(gffline)
        
        transcript_set.add(gffline['transcript_id'])
        
        #temp_gff_list.append(gffline)    
    
    ## DEBUG
    #~ print "processed %s lines in %s " % (i, gff_fn)
    #~ print "obtained %s transcript names" % (len(transcript_set))
    #~ print "values/exon types counts: "
    #~ print gff_values_counts
    
    exon_numbers       = {}
    plus_strand_trans  = 0
    minus_strand_trans = 0
    for transcript_id in transcript_dict.keys():
        current_record = transcript_dict[transcript_id]
        num_of_exons = len(current_record)
        #print num_of_exons
        if num_of_exons not in exon_numbers.keys():
            exon_numbers[num_of_exons] = 1
        else:
            exon_numbers[num_of_exons] += 1
        
        if  current_record[0]['strand'] == '+':
            plus_strand_trans += 1
        
        elif current_record[0]['strand'] == '-':
            minus_strand_trans += 1
        else:
            print "bad record, no strand: ", current_record
    ## debug
    ## print "plus_strand_trans: %s " % (plus_strand_trans)
    ## print "minus_strand_trans: %s " % (minus_strand_trans)
    
    
    ## For future gene selection & DEBUG
    ##for key in exon_numbers.keys():
    ##    print "%s exons: %s transcripts" % (  key, exon_numbers[key] )
    
    
    return transcript_dict
    
    
    
    



def gfflines(filename):
    """Open an optionally gzipped GTF file and generate a dict for each line.
    """
    #fn_open = gzip.open if filename.endswith('.gz') else open

    with open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield gffparse(line)
            
                

def gffparse(line):
    """Parse a single GTF line and return a dict.
    """
    GTF_HEADER = cfg['GTF_HEADER']
    result = {}

    fields = line.rstrip().split('\t')
    if len(fields) !=9:
        print line,
        print fields
    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])
    return result
    
def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value
    
    



##debug
## print transcript_dict

genome_fas = Fasta( cfg['genome_fasta'], sequence_always_upper=True)
genome_fas_dict = {}
for key in genome_fas.keys():
    genome_fas_dict[key] = len(genome_fas[key][:])

#print "got here"
#print  '10274.m000709', genome_fas_dict['10274.m000709']
#extract sequences

def extract_splices(transcript_dict, genome_fas, genome_fas_dict, splice_type):

    def plus_strand_feat(chrom, position, offset):
        chrom_size = genome_fas_dict[chrom]
        #print chrom, chrom_size, position, offset
        prefix_len  = 0
        postfix_len = 0
        prefix_str  = ""
        postfix_str = ""
        
        if (position  - offset)  < 1:
            prefix_len = abs(position - offset) + 1
            prefix_str = "N" * prefix_len
            temp_seq = "%s%s" % (prefix_str, genome_fas[chrom][:position + offset +1])
        elif position+offset > chrom_size:
            postfix_len = (position + offset) - chrom_size
            postfix_str = "N" * postfix_len
            temp_seq = "%s%s" % (genome_fas[chrom][position - offset:chrom_size], postfix_str)
        else:
            temp_seq = genome_fas[chrom][position - offset-1:position + offset-1]
        return temp_seq, len(temp_seq)
        
    def minus_strand_feat(chrom, position, offset):
        chrom_size = genome_fas_dict[chrom]
        #print chrom, chrom_size, position, offset
        prefix_len  = 0
        postfix_len = 0
        prefix_str  = ""
        postfix_str = ""
        
        if (position  - offset)  < 1:
            prefix_len = abs(position-offset) + 1 
            prefix_str = "N" * prefix_len
            temp_seq = prefix_str + genome_fas[chrom][:position+offset+1].complement
        elif position+offset > chrom_size:
            postfix_len = (position+offset) - chrom_size
            postfix_str = "N" * postfix_len
            temp_seq = genome_fas[chrom][position-offset:chrom_size].complement + postfix_str
        else:
            temp_seq = genome_fas[chrom][position-offset:position+offset].complement
        temp_seq = temp_seq[::-1] #reversing
        return temp_seq, len(temp_seq)
    
    
    #counter = 0
    def get_acceptors(transcript_list):
        if transcript_list == None:
            transcript_list = transcript_dict.keys()
        else:
            pass
            ## DEBUG
            ##print transcript_list
     
        for transcript_id in transcript_list:
            current_record = transcript_dict[transcript_id]
            num_of_exons = len(current_record)
            #print "###", transcript_id, num_of_exons
            #print transcript_id, num_of_exons
            
            if num_of_exons > 1:
                strand = current_record[0]['strand']
                if strand == '+':
                    for dict_exon_num in range(1, num_of_exons):
                        chrom    = current_record[dict_exon_num]['seqname']
                        position = current_record[dict_exon_num]['start']
                        seq, seq_len = plus_strand_feat(chrom, position, offset)
                        real_exon_num = dict_exon_num + 1
                        print "%s_Fex%s_%s_%s_%s__%s\t%s" % (transcript_id, real_exon_num, label, seq_len, seq[offset-2:offset], position, seq)
                elif strand == '-':
                    for dict_exon_num in range(0, num_of_exons-1):
                        chrom    = current_record[dict_exon_num]['seqname']
                        position = current_record[dict_exon_num]['end']
                        seq, seq_len = minus_strand_feat(chrom, position, offset)
                        real_exon_num = num_of_exons - dict_exon_num
                        print "%s_Rex%s_%s_%s_%s__%s\t%s" % (transcript_id, real_exon_num, label, seq_len, seq[offset-2:offset], position, seq)
                else:
                    print "XXX strand error: no strand info for transcript_id"
                    print "XXX", current_record
    def get_donors(transcript_list):
        if transcript_list == None:
            transcript_list = transcript_dict.keys()
     
        for transcript_id in transcript_list:
            current_record = transcript_dict[transcript_id]
            num_of_exons = len(current_record)
            print "###", transcript_id, num_of_exons
            
            if num_of_exons > 1:
                strand = current_record[0]['strand']
                if strand == '+':
                    for dict_exon_num in range(0, num_of_exons-1):
                        chrom    = current_record[dict_exon_num]['seqname']
                        position = current_record[dict_exon_num]['end']
                        seq, seq_len = plus_strand_feat(chrom, position, offset)
                        real_exon_num = dict_exon_num + 1
                        print "%s_Fex%s_%s_%s_%s__%s\t%s" % (transcript_id, real_exon_num, label, seq_len, seq[offset+1:offset+3],position, seq)
                elif strand == '-':
                    for dict_exon_num in range(1, num_of_exons):
                        chrom    = current_record[dict_exon_num]['seqname']
                        position = current_record[dict_exon_num]['start']
                        seq, seq_len = minus_strand_feat(chrom, position, offset)
                        real_exon_num = num_of_exons - dict_exon_num
                        print "%s_Rex%s_%s_%s_%s__%s\t%s" % (transcript_id, real_exon_num, label, seq_len, seq[offset+1:offset+3], position, seq)
                else:
                    print "XXX strand error: no strand info for transcript_id"
                    print "XXX", current_record            
        #the AG from

    def get_ATGs(transcript_list):
        if transcript_list == None:
            transcript_list = transcript_dict.keys()
     
        for transcript_id in transcript_list:
            current_record = transcript_dict[transcript_id]
            chrom    = current_record[0]['seqname']
            num_of_exons = len(current_record)
            strand = current_record[0]['strand']
            if strand == "+":
                position = current_record[0]['start']
                
                seq, seq_len = plus_strand_feat(chrom, position, offset)
                print "%s_Fstart_%s__%s\t%s" % (transcript_id, seq[offset:offset+3], position, seq)
            if strand == "-":
                position = current_record[-1]['end']
                
                seq, seq_len = minus_strand_feat(chrom, position, offset)
                print "%s_Rstart_%s__%s\t%s" % (transcript_id, seq[offset:offset+3], position, seq)
    
    def get_stops(transcript_list):
        if transcript_list == None:
            transcript_list = transcript_dict.keys()
     
        for transcript_id in transcript_list:
            current_record = transcript_dict[transcript_id]
            chrom    = current_record[0]['seqname']
            num_of_exons = len(current_record)
            strand = current_record[0]['strand']
            if strand == "+":
                position = current_record[-1]['end'] -2
                
                seq, seq_len = plus_strand_feat(chrom, position, offset)
                print "%s_Fstop_%s__%s\t%s" % (transcript_id, seq[offset:offset+3], position, seq)
            if strand == "-":
                position = current_record[0]['start'] +2
                
                seq, seq_len = minus_strand_feat(chrom, position, offset)
                print "%s_Rstop_%s__%s\t%s" % (transcript_id, seq[offset:offset+3], position, seq)
    
    #executing the functions
    offset   = 30
    position = 0
    label = splice_type[:3]
    transcript_list = transcript_dict.keys()
    
    if  splice_type == "acceptors":
        #AG types, exon numbering  from 0
        #assuming position sorted not nested genes
                #F_strand_exons_num_range = ('start', 1, "last")
        #R_strand_exons_num_range = ('end',   0, "second_last")
        get_acceptors(transcript_list)
        
    elif  splice_type == "donors":
        get_donors(transcript_list)
    elif  splice_type == "ATG":
        get_ATGs(transcript_list) 
    elif  splice_type == "stop":
        get_stops(transcript_list)   
    else:
        print "ERROR: not implemented yet"


def extract_cds_seqs(*transcript_list):
    transcript_list = transcript_dict.keys()
    #debug
    #transcript_list = ['model.1879.m002022', 'foobar']
    #pattern for longest ORF in transcript
    orf_finder = re.compile(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)')
    
    def get_exon_seq(exon_record, strand):            
        current_exon = exon_record
        chrom    = current_exon['seqname']
        start    = current_exon['start']
        end      = current_exon['end']
        if strand == "+":                
            exon_seq = genome_fas[chrom][start-1:end]
        if strand == "-":
            exon_seq = genome_fas[chrom][start-1:end] #.complement
            #exon_seq = exon_seq[::-1]
        #print "RRR", len(exon_seq)
        #print "RRR\t%s" % (exon_seq)
        return exon_seq
        
    
    
    for transcript_id in transcript_list:
        current_record = transcript_dict[transcript_id]
        num_of_exons = len(current_record)
        strand = current_record[0]['strand']
        
        transcript_seq = ""
        cds_seq        = "" 
        if num_of_exons == 1:
            current_exon = current_record[0]
            transcript_seq = "%s" % get_exon_seq(current_exon, strand)
        else:
            for exon in current_record:
                #print "$$$$", exon
                transcript_seq += "%s" % get_exon_seq(exon, strand)
        if strand == '-':
            tmp_seq =  Seq(transcript_seq)
            transcript_seq = "%s" %(tmp_seq.reverse_complement())

        ## DEBUG
        #~ print ">%s_from_%s_exons_%s_trans" % (transcript_id, num_of_exons, strand)
        #~ print "%s" % transcript_seq
        try:        
            cds_seq = max(orf_finder.findall(transcript_seq), key = len)
            print ">%s_from_%s_exons_%s_CDS" % (transcript_id, num_of_exons, strand)
            print cds_seq
        except:
            print "error: transcript_id", transcript_id

def extract_introns(*transcript_list):
    transcript_list   = transcript_dict.keys()
    intron_sizes_list = []
    
    #~ def get_exon_seq(exon_record, strand):            
        #~ current_exon = exon_record
        #~ chrom    = current_exon['seqname']
        #~ start    = current_exon['start']
        #~ end      = current_exon['end']
        #~ if strand == "+":                
            #~ exon_seq = genome_fas[chrom][start-1:end]
        #~ if strand == "-":
            #~ exon_seq = genome_fas[chrom][start-1:end] #.complement
            #~ #exon_seq = exon_seq[::-1]
        #~ #print "RRR", len(exon_seq)
        #~ #print "RRR\t%s" % (exon_seq)
        #~ return exon_seq
        
    
    
    for transcript_id in transcript_list:
        current_record = transcript_dict[transcript_id]
        num_of_exons = len(current_record)
        if num_of_exons != 1:
            introns_positions_list = []
            chrom = current_record[0]['seqname']
            strand = current_record[0]['strand']
            introns_seq = ""
            for dict_exon_num in range(0, num_of_exons):
                if dict_exon_num == 0:
                    intron_start = current_record[dict_exon_num]['end'] + 1
                elif dict_exon_num == num_of_exons -1:
                    intron_end = current_record[dict_exon_num]['start'] - 1
                    introns_positions_list.append([intron_start, intron_end])
                else:
                    intron_end = current_record[dict_exon_num]['start'] - 1
                    introns_positions_list.append([intron_start, intron_end])
                    intron_start = current_record[dict_exon_num]['end'] + 1
            ##print "INTR:",     introns_positions_list
            for intron_pos in introns_positions_list:
                start, end = intron_pos
                seq = genome_fas[chrom][start-1:end]
                intron_len = len(seq)
                print "INTR_LEN:", intron_len
                intron_sizes_list.append(intron_len)
                ##DEBUG
                introns_seq = '%s\n%s' % (introns_seq, seq)
                ##introns_seq = '%s%s' % (introns_seq, seq)
            print ">%s_introns_%s__%s_%s_%s" % (transcript_id, num_of_exons, current_record[0]['end'] +1, current_record[-1]['start'] -1, strand)
            print introns_seq
                #print seq
        else:
            pass
    
    intron_sizes_list.sort()
    print "INTRON_MIN_MAX",    intron_sizes_list[0], intron_sizes_list[-1] 
        #~ for exon in current_record:
                #~ #print "$$$$", exon
                #~ introns_seq += "%s" % get_intron_seq(exon, strand)
        #~ if strand == '-':
            #~ tmp_seq =  Seq(transcript_seq)
            #~ transcript_seq = "%s" %(tmp_seq.reverse_complement())

        #~ ## DEBUG
        #~ print ">%s_from_%s_exons_%s_trans" % (transcript_id, num_of_exons, strand)
        #~ print "%s" % transcript_seq
        #~ try:        
            #~ cds_seq = max(orf_finder.findall(transcript_seq), key = len)
            #~ print ">%s_from_%s_exons_%s_CDS" % (transcript_id, num_of_exons, strand)
            #~ print cds_seq
        #~ except:
            #~ print "error: transcript_id", transcript_id


if __name__ == "__main__":
    gff_fn = cfg['gff']['genes_fn']
    transcript_dict = gff_reader(gff_fn)
    #print "transcript_dict size:", len(transcript_dict)
    
    extract_splices(transcript_dict, genome_fas, genome_fas_dict, "acceptors")
    extract_splices(transcript_dict, genome_fas, genome_fas_dict, "donors")
    extract_cds_seqs()
    extract_splices(transcript_dict, genome_fas, genome_fas_dict, "ATG")
    extract_splices(transcript_dict, genome_fas, genome_fas_dict, "stop")
    extract_introns()
