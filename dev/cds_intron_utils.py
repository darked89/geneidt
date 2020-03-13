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
        exon_seq = genome_fas[chrom][start-1:end]
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
            print(">%s_from_%s_exons_%s_CDS" % (transcript_id, num_of_exons, strand))
            print( cds_seq)
        except:
            print("error: transcript_id", transcript_id)
            
            

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
            for index in range(0, len(introns_positions_list)):
                intron_pos = introns_positions_list[index]
                start, end = intron_pos
                intron_seq_name = "%s_intron_%s__%s_%s" % (transcript_id, index+1, start, end)
                intron_seq = genome_fas[chrom][start-1:end]
                if strand == '-':
                    tmp_seq =  Seq('%s' % intron_seq)
                    intron_seq = "%s" % (tmp_seq.reverse_complement())
                print("%s\t%s" % (intron_seq_name, intron_seq))
                intron_len = len(intron_seq)
                #print "INTR_LEN:", intron_len
                intron_sizes_list.append(intron_len)
                ##DEBUG
                
                ##introns_seq = '%s\n%s' % (introns_seq, seq)
                ##introns_seq = '%s%s' % (introns_seq, seq)
            ##print ">%s_introns_%s__%s_%s_%s" % (transcript_id, num_of_exons, current_record[0]['end'] +1, current_record[-1]['start'] -1, strand)
            ##print introns_seq
                #print seq
        else:
            pass
    
    intron_sizes_list.sort()
    print("INTRON_MIN_MAX",    intron_sizes_list[0], intron_sizes_list[-1])
            
