#!/usr/bin/env python3

import sys, os

from random import randint

from pyfaidx import Fasta

kmer_size = 60

fasta_in = sys.argv[1]
tbl_out  = sys.argv[2]
out_fas_dir  = sys.argv[3]

#print "Hallo from the py script"

genome_fas = Fasta( fasta_in, sequence_always_upper=True)

out_fh = open(tbl_out, "w")

for contig_name in genome_fas.keys():
    out_seq = genome_fas[contig_name][:]
    changed_seq_name = contig_name.replace(".","")
    out_str = "%s\t%s\n" % (changed_seq_name, out_seq)
    out_fh.write(out_str)
    out_fas_fn  = "%s/%s" % (out_fas_dir, changed_seq_name)
    out_fas_fh   = open(out_fas_fn, "w")
    my_string = ">%s\n" % (contig_name)
    out_fas_fh.write(my_string)
    my_string = "%s\n" % (out_seq)
    out_fas_fh.write(my_string)
    out_fas_fh.close()

    out_len_fname = out_fas_fn + "_len"
    #print out_len_fname
    out_len_fh = open(out_len_fname, "w")
    my_string = "%s\t%s\n" % (changed_seq_name, len( out_seq))
    out_len_fh.write(my_string);
    out_len_fh.close()

out_fh.close()
