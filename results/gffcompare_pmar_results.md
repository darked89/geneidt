
```
# whole genome prediction
~/proj/geneidc/bin/geneid -G -3 -P pmar01.geneid.optimized.param ./../genomes_4geneidc_t/p.mari/Perkinsus_marinus_atcc_50983_gca_000006405.JCVI_PMG_1.0.dna.toplevel.fa > geneid_test_genome_pmar.gff

# hack to use geneid_1.3.15dev gff3 output with gffcompare

sed -i 's/CDS/exon/g' geneid_test_genome_pmar.gff

```

# gffcompare v0.11.6 | Command line was:
#/home/darked/soft/bin/gffcompare -V -r ../../genomes_4geneidc_t/p.mari/Perkinsus_marinus_atcc_50983_gca_000006405.JCVI_PMG_1.0.46.gtf -o genome_pmar geneid_test_genome_pmar.gff
#

#= Summary for dataset: geneid_test_genome_pmar.gff 
#     Query mRNAs :   31773 in   31773 loci  (17895 multi-exon transcripts)
#            (0 multi-transcript loci, ~1.0 transcripts per locus)
# Reference mRNAs :   37681 in   37497 loci  (17104 multi-exon)
# Super-loci w/ reference transcripts:    17955
#-----------------| Sensitivity | Precision  |
        Base level:    68.4     |    73.5    |
        Exon level:    40.6     |    47.8    |
      Intron level:    45.6     |    53.6    |
Intron chain level:     9.1     |     8.7    |
  Transcript level:     8.5     |    10.1    |
       Locus level:     8.5     |    10.1    |

     Matching intron chains:    1558
       Matching transcripts:    3202
              Matching loci:    3202

          Missed exons:   53199/145812	( 36.5%)
           Novel exons:   34267/123772	( 27.7%)
        Missed introns:   24576/108127	( 22.7%)
         Novel introns:   20438/91999	( 22.2%)
           Missed loci:   15946/37497	( 42.5%)
            Novel loci:   13083/31773	( 41.2%)

 Total union super-loci across all input datasets: 31038 
31773 out of 31773 consensus transcripts written in debug_pmar.annotated.gtf (0 discarded as redundant)
