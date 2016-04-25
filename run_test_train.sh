#!/bin/sh

echo "running: run_test_train.sh with pmarinus4training data";
./geneid_trainer.pl -species pmar01 -gff test_data/pmarinus4training.gff -fastas test_data/pmarinus4training.fa \
-sout statsfile_out.pmar01
