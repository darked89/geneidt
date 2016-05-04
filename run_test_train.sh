#!/bin/sh

echo "running: run_test_train.sh with pmarinus4training data";

#./gff2geneid_gff.pl -gff test_data/pmarinus4training.gff


./geneid_trainer.pl -species pmar01 -gff test_data/pmarinus4training.gff -fasta test_data/pmarinus4training.fa \
-sout statsfile_out.pmar01 -fix_fasta 0
