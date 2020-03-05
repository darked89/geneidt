
```
# whole genome prediction
~/proj/geneidc/bin/geneid -G -3 -P pmar01.geneid.optimized.param ./../genomes_4geneidc_t/p.mari/Perkinsus_marinus_atcc_50983_gca_000006405.JCVI_PMG_1.0.dna.toplevel.fa > geneid_test_genome_pmar.gff

# hack to use geneid_1.3.15dev gff3 output with gffcompare

sed -i 's/CDS/exon/g'

```

# gffcompare v0.11.6 | Command line was:
#/home/darked/soft/bin/gffcompare -V -r ../../genomes_4geneidc_t/p.mari/Perkinsus_marinus_atcc_50983_gca_000006405.JCVI_PMG_1.0.46.gtf -o genome_pmar geneid_test_genome_pmar.gff
#

#= Summary for dataset: geneid_test_genome_pmar.gff 
#     Query mRNAs :   31987 in   31987 loci  (17890 multi-exon transcripts)
#            (0 multi-transcript loci, ~1.0 transcripts per locus)
# Reference mRNAs :   37681 in   37497 loci  (17104 multi-exon)
# Super-loci w/ reference transcripts:    17978
#-----------------| Sensitivity | Precision  |
        Base level:    67.9     |    72.7    |
        Exon level:    39.2     |    46.8    |
      Intron level:    43.5     |    52.1    |
Intron chain level:     8.7     |     8.3    |
  Transcript level:     8.4     |     9.9    |
       Locus level:     8.4     |     9.9    |

     Matching intron chains:    1481
       Matching transcripts:    3157
              Matching loci:    3157

          Missed exons:   54747/145812	( 37.5%)
           Novel exons:   34534/122292	( 28.2%)
        Missed introns:   23820/108127	( 22.0%)
         Novel introns:   20512/90305	( 22.7%)
           Missed loci:   15644/37497	( 41.7%)
            Novel loci:   13388/31987	( 41.9%)

 Total union super-loci across all input datasets: 31366 
31987 out of 31987 consensus transcripts written in 1092efo3fixed_pmar.annotated.gtf (0 discarded as redundant)
