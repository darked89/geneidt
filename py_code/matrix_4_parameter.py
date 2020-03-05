#!/usr/bin/env python3
"""
#!/usr/bin/gawk -f
#USAGE= $BIN/preparetrimatrixstart4parameter.awk 4 5 6 7 $DATA/set2/sites/starts.geneid.dimatrix
BEGIN {
    st = ARGV[1];
    nd = ARGV[2];
    rd = ARGV[3];
    th = ARGV[4];   
    ARGV[1] = ARGV[2] = ARGV[3] = ARGV[4] = "";
}
{
  if (($1 == st && $2 !~ /A$/) || 
      ($1 == nd && $2 != "AT") || 
      ($1 == rd && $2 != "TG") || 
      ($1 == th && $2 !~ /^G/))
      {
         print $1, $2, -9999;
    } else 
    if (($1 == st && $2 ~  /A$/) ||\
        ($1 == nd && $2 == "AT") || 
        ($1 == rd && $2 == "TG"))
    print $1, $2, 0;
  else
print;
}
"""

import sys
from itertools import product

input_fn  = sys.argv[1]
site_type = sys.argv[2]
positions = sys.argv[3:6]

positions = [int(x) for x in positions]

## DEBUG
#~ print input_fn, site_type, positions

sites_dict = {'donor' : 'GT',
              'acceptor' : 'AG',
              ##'start' : 'ATG',
}

def read_3_cols(input_fn):
    tmp_logratio_dict = {}
    result_freq_dict = {}
    for line in open(input_fn).readlines():
        sl = line.split()
        position  = int(sl[0])
        kmer      = sl[1]
        logratio = float(sl[2])
        if position not in tmp_logratio_dict.keys():
            tmp_logratio_dict[position] = {}
        tmp_logratio_dict[position][kmer] = logratio
    ##result_freq_dict 
    return tmp_logratio_dict
    
four_bases = 'ACGT'
kmer_plus_list2 = []
for element in product(four_bases, repeat=2):
    my_string = ''.join(x for x in element)
    #print my_string
    kmer_plus_list2.append(my_string)
    

site_logratios_dict = read_3_cols(input_fn)
motif_bases = sites_dict[site_type]

for position in positions:
    #position = int(position)
    for kmer in site_logratios_dict[position].keys():
        list_index = positions.index(position)
        if list_index == 0: 
            if kmer[-1] != motif_bases[0]:
                site_logratios_dict[position][kmer] = -9999
            else:
                pass
        elif list_index == 1:
            if kmer != motif_bases[:2]:
                site_logratios_dict[position][kmer] = -9999
            else:
                site_logratios_dict[position][kmer] = 0
        elif list_index == 2:
                if kmer[0] != motif_bases[1]:
                    site_logratios_dict[position][kmer] = -9999

for matrix_position in site_logratios_dict.keys():
    for kmer in kmer_plus_list2:
        logratio = site_logratios_dict[matrix_position][kmer]
        out_string = '%s %s %s' % (matrix_position, kmer, logratio)
        print( out_string)
