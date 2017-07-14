

def extract_sites(transcript_list, site_type, transcript_dict, genome_fas, genome_fas_dict):  
    ##transcript_list = transcript_dict.keys()
    sites_dic = {}
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
        for transcript_id in transcript_list:
            current_record = transcript_dict[transcript_id]
            num_of_exons = len(current_record)
            #print "###", transcript_id, num_of_exons
            #print transcript_id, num_of_exons
            
            if num_of_exons > 1:
                strand = current_record[0]['strand']
                chrom  = current_record[0]['seqname']
                if strand == '+':
                    for dict_exon_num in range(1, num_of_exons):
                        position = current_record[dict_exon_num]['start']
                        seq, seq_len = plus_strand_feat(chrom, position, offset)
                        real_exon_num = dict_exon_num + 1
                        key = "%s_Fex%s_%s_%s_%s__%s" % (transcript_id, real_exon_num, label, seq_len, seq[offset-2:offset], position)
                        sites_dic[key] = '%s' % (seq) 
                        #print("%s_Fex%s_%s_%s_%s__%s\t%s" % (transcript_id, real_exon_num, label, seq_len, seq[offset-2:offset], position, seq))
                elif strand == '-':
                    for dict_exon_num in range(0, num_of_exons-1):
                        position = current_record[dict_exon_num]['end']
                        seq, seq_len = minus_strand_feat(chrom, position, offset)
                        real_exon_num = num_of_exons - dict_exon_num
                        key = "%s_Rex%s_%s_%s_%s__%s" % (transcript_id, 
                        real_exon_num, label, seq_len, seq[offset-2:offset], position)
         
                        sites_dic[key] = '%s' % (seq) 
                        #print("%s_Rex%s_%s_%s_%s__%s\t%s" % (transcript_id, real_exon_num, label, seq_len, seq[offset-2:offset], position, seq))
                else:
                    print( "XXX strand error: no strand info for transcript_id")
                    print("XXX", current_record)
    def get_donors(transcript_list):
     
        for transcript_id in transcript_list:
            current_record = transcript_dict[transcript_id]
            num_of_exons = len(current_record)
            print("###", transcript_id, num_of_exons)
            
            if num_of_exons > 1:
                strand = current_record[0]['strand']
                chrom  = current_record[0]['seqname']
                if strand == '+':
                    for dict_exon_num in range(0, num_of_exons-1):
                        
                        position = current_record[dict_exon_num]['end']
                        seq, seq_len = plus_strand_feat(chrom, position, offset)
                        real_exon_num = dict_exon_num + 1
                        key = "%s_Fex%s_%s_%s_%s__%s" % (transcript_id, real_exon_num, label, seq_len, seq[offset+1:offset+3],position)
                        sites_dic[key] = '%s' % (seq) 
                        #print("%s_Fex%s_%s_%s_%s__%s\t%s" % (transcript_id, real_exon_num, label, seq_len, seq[offset+1:offset+3],position, seq))
                elif strand == '-':
                    for dict_exon_num in range(1, num_of_exons):
                        
                        position = current_record[dict_exon_num]['start']
                        seq, seq_len = minus_strand_feat(chrom, position, offset)
                        real_exon_num = num_of_exons - dict_exon_num
                        key = "%s_Rex%s_%s_%s_%s__%s" % (transcript_id, real_exon_num, label, seq_len, seq[offset+1:offset+3],position)
                        sites_dic[key] = '%s' % (seq) 
                        #print("%s_Rex%s_%s_%s_%s__%s\t%s" % (transcript_id, real_exon_num, label, seq_len, seq[offset+1:offset+3], position, seq))
                else:
                    print("XXX strand error: no strand info for transcript_id")
                    print("XXX", current_record)            
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
                key = "%s_Fstart_%s__%s" % (transcript_id, seq[offset:offset+3], position)
                sites_dic[key] = '%s' % (seq)
                #print("%s_Fstart_%s__%s\t%s" % (transcript_id, seq[offset:offset+3], position, seq))
            if strand == "-":
                position = current_record[-1]['end']
                
                seq, seq_len = minus_strand_feat(chrom, position, offset)
                key = "%s_Rstart_%s__%s" % (transcript_id, seq[offset:offset+3], position)
                sites_dic[key] = '%s' % (seq)
                ##print("%s_Rstart_%s__%s\t%s" % (transcript_id, seq[offset:offset+3], position, seq))
    
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
                key = "%s_Fstop_%s__%s" % (transcript_id, seq[offset:offset+3], position)
                sites_dic[key] = '%s' % (seq)
                ##print("%s_Fstop_%s__%s\t%s" % (transcript_id, seq[offset:offset+3], position, seq))
            if strand == "-":
                position = current_record[0]['start'] +2
                
                seq, seq_len = minus_strand_feat(chrom, position, offset)
                key = "%s_Rstop_%s__%s" % (transcript_id, seq[offset:offset+3], position)
                sites_dic[key] = '%s' % (seq)
                
                ##print("%s_Rstop_%s__%s\t%s" % (transcript_id, seq[offset:offset+3], position, seq))
    
    #executing the functions
    offset   = 30
    position = 0
    label = site_type[:3]
    revforward_labels = {'+': 'F', '-':'R'}
    
    
    if  site_type == "acceptors":
        get_acceptors(transcript_list)
    elif  site_type == "donors":
        get_donors(transcript_list)
    elif  site_type == "ATG":
        get_ATGs(transcript_list) 
    elif  site_type == "stop":
        get_stops(transcript_list)   
    else:
        print("ERROR: not implemented yet")
    return sites_dic


def base_freq_count(seq_dict):
    nucleotides_lst  = ('A', 'C', 'G', 'T', 'N')
    nucleotides_clean  = ('A', 'C', 'G', 'T' ) 
    counter = 0
    while counter < 1:
        for key in seq_dict.keys():
            seq_len = len(seq_dict[key])
            counter +=1

    print(seq_len)
    freq_dict = {}
    #print new_seq_len
    for base_pos in range(0, seq_len):
        freq_dict[base_pos] = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'ACGTsum': 0} 
        for key in seq_dict.keys():
             current_base = seq_dict[key][base_pos]
             assert current_base in nucleotides_lst
             if current_base != 'N':
                freq_dict[base_pos][current_base] += 1
                freq_dict[base_pos]['ACGTsum']    += 1

    freq_output_list = []
    for position in freq_dict.keys():
		#tmp_position_dict = {}
		#freq_output_dict[position] = {}
        for base in nucleotides_clean:
            base_count       = freq_dict[position][base]
            total_base_count = freq_dict[position]['ACGTsum']
            base_frequency   = base_count / total_base_count
            freq_output_list.append([position+1, base, base_frequency, base_count])
            #OK working line below
            #print("%s\t%s\t%s\t%.4f" % (position+1, base, base_frequency, base_count))
    return  freq_output_list           
