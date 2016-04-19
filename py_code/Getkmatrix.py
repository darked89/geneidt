#!/usr/bin/env python2

'''
translating:

# Getkmatrix.awk
# gpf, imim, January 2000
# usage: gawk -f Getkmatrix.awk order number_of_nucleotides_per_sequence cds_seqfile (.tbl)

BEGIN {
  PCOUNT  = 0.25;
  k       = ARGV[1];          #order of the markov chain
  num     = ARGV[2];
  ARGV[1] = ARGV[2] = "";

  # set pseudocounts 
  alpha[1] = "A";
  alpha[2] = "C";
  alpha[3] = "G";
  alpha[4] = "T";

  # for all k_tuples 
  for (i=1; i<=4; i++) 
     nx(1, k, alpha[i],  N0, 4*PCOUNT, num-k);

  # for all k+1_tuples
  for (i=1; i<=4; i++) 
     nx(1, k+1, alpha[i], N,   PCOUNT, num-k);

#  sizek1=4^(k+1);

}

{
  sequence = toupper($2);
  if (sequence !~ /[^ACGT]/) { # consider only standard acgt

    lseq = length(sequence); 
    L   += (lseq-k);

    for (i=1; i<=lseq-k; i++) {
      ktuple   = substr(sequence, i, k);
      ktuple1  = substr(sequence, i, k+1);
      
      N0[i, ktuple]++; # ktuple frequence at position i
      N[i, ktuple1]++; # ktuple1 frequence at position i
      
    }
    total_seq++
  }
}
END {
  # get number of k-tuples observed in each frame
  L0 = L+(3*sizek);

  for (t in N) {        # transition probabilities
    split (t, x, SUBSEP);
    pos = x [1];
    tk  = substr(x[2], 1, k);
    tk1 = x[2];
    print x[1], x[2], N[t]/N0[pos,tk];
  }
}

function nx(l, len, s, Mat, p, num, i, pos) { 
#nx(1, k, alpha[i],  N0, 4*PCOUNT, num-k);
  if (l == len) {
    for (pos=1; pos<=num; pos++)
      Mat[pos,s] = p;
  }
  else {
     l++;
     for (i=1; i<=4; i++) 
       nx(l, len, s alpha[i], Mat, p, num);
   }
}

'''

import sys
from itertools import product


PCOUNT  = 0.25
DNAbases = ['A', 'C', 'G', 'T']

kmer_size  = int(sys.argv[1])          #order of the markov chain
seq_len_in = int(sys.argv[2])
in_fn      = sys.argv[3]

kmer_plus_size = kmer_size + 1 # we are dealing with dimers here so far
kmer_plus_list = []
for base_A in DNAbases:
	for base_B in DNAbases:
		tmp_dimer =  base_A + base_B
		kmer_plus_list.append(tmp_dimer)

print kmer_plus_list

four_bases = 'ACGT'
kmer_plus_list2 = []
for element in product(four_bases, repeat=2):
	my_string = ''.join(x for x in element)
	#print my_string
	kmer_plus_list2.append(my_string)


def nx(init_kmer_size, kmer_size, base, frequence_dict,   PCOUNT, last_kmer_pos):
    if init_kmer_size == kmer_size:
		for i in range(1, len(frequence_dict.keys()):
			
	else:
			
    
    #x(init_kmer_size, kmer_size, base, matrix_A,   PCOUNT, last_kmer_pos)
'''
for base in DNAbases:
	nx(1, kmer_size,    base, 'N0', 4*PCOUNT, seq_len-kmer_size)
    nx(1, kmer_size +1, base, 'N',   PCOUNT,  seq_len-kmer_size ) #error here??
    
    
'''

ktuple_dict      = {}
ktuple_plus_dict = {}

for line in open(in_fn).readlines():
	sl = line.split()
	sequence = sl[1]
	discard = 0
	for base in sequence:
		if base not in DNAbases:
			discard = 1
			break
	if discard != 1:
		seq_len_real = len(sequence)
		for i in range(0, seq_len_real - kmer_size + 1):
			tmp_kmer = sequence[i:i+kmer_size]
			position = i+1
			print position, tmp_kmer
			if position not in ktuple_dict.keys():
				ktuple_dict[position] = {}
				for base in DNAbases:
					ktuple_dict[position][base] = 0
					#{'A':0, 'C':0, 'G':0, 'T':0}
			ktuple_dict[position][tmp_kmer] += 1 #increment count
		
		kmer_plus_size = kmer_size + 1
		for i in range(0, seq_len_real - kmer_plus_size + 1):
			tmp_kmer = sequence[i:i+kmer_plus_size]
			position = i+1
			if position not in ktuple_plus_dict.keys():
				ktuple_plus_dict[position] = {}
				for list_kmer in kmer_plus_list2:
					ktuple_plus_dict[position][list_kmer] = 0
			ktuple_plus_dict[position][tmp_kmer]  += 1 #increment count
	else:
		pass 
		#one can maybe fix this, because we just need to discard N/etc. containing tuples, not the whole seq...
		
for position in ktuple_dict.keys():
	print position, ktuple_dict[position]
	

for position in ktuple_plus_dict.keys():
	print position, ktuple_plus_dict[position]
