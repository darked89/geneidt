
def gff_reader(gff_fn, cfg):
    #~ GTF_HEADER = cfg['GTF_HEADER']
    #~ #temp_gff_list     = []
    #~ transcript_set = set()
    #~ transcript_dict = {}
    #~ gff_values_counts = {}
    #~ int_values =  cfg['gff']['int_values']
    #~ col3_vals = cfg['gff']['col3_vals']

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
        
        result = {}   
        fields = line.rstrip().split('\t')
        if len(fields) !=9:
            print(line,)
            print(fields)
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

    #TEST
    GTF_HEADER = cfg['GTF_HEADER']
    #temp_gff_list     = []
    transcript_set = set()
    transcript_dict = {}
    gff_values_counts = {}
    int_values =  cfg['gff']['int_values']
    col3_vals = cfg['gff']['col3_vals']
    #execution starts here
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
            print("bad record, no strand: %" % current_record)
    ## debug
    ## print "plus_strand_trans: %s " % (plus_strand_trans)
    ## print "minus_strand_trans: %s " % (minus_strand_trans)
    
    
    ## For future gene selection & DEBUG
    ##for key in exon_numbers.keys():
    ##    print "%s exons: %s transcripts" % (  key, exon_numbers[key] )
    
    
    return transcript_dict
    
    
    
    



            
                


