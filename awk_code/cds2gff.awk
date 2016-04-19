#!/usr/bin/gawk -f 

BEGIN{
  # default source name
  source="ANNOTATED";

  # data structure
  locus  = 1;
  seqlen = 2;
  strand = 3;
  fcds   = 4;  # first cds position

  # gff features
  FIRST    = "First";
  INTERNAL = "Internal";
  TERMINAL = "Terminal";
  SINGLE   = "Single";

  # gff defaults
  SCORE = 1;
  GROUP = 1;

  # gff formatted output
  gff_format="%s\t%s\t%s\t%d\t%d\t%f\t%s\t%d\t%s\n";
}
{
  frame=0;
  group[$locus]++;

  # check if single
  if (NF==fcds+1) {
      printf gff_format, $locus, source, SINGLE, $fcds, $(fcds+1), SCORE, $strand, frame, $locus"_"group[$locus];  
      next;
  }

  # assumes first first exon, last last exon, others internal

  # print first exon
  remainder=($(fcds+1)-($fcds+frame)+1)%3;
  nextframe=(3-($(fcds+1)-($fcds+frame)+1)%3)%3;
  
  if ($strand == "+") 
    printf gff_format, $locus, source, FIRST, $fcds, $(fcds+1), SCORE, $strand, frame, $locus"_"group[$locus];  
  else
    printf gff_format, $locus, source, TERMINAL, $fcds, $(fcds+1), SCORE, $strand, remainder, $locus"_"group[$locus];  

  # print internal exons
  for (i=fcds+2;i<NF-2;i+=2) {
    frame=nextframe;
    remainder=($(i+1)-($i+frame)+1)%3;
    nextframe=(3-($(i+1)-($i+frame)+1)%3)%3;
  if ($strand == "+")
    printf gff_format, $locus, source, INTERNAL, $i, $(i+1), SCORE, $strand, frame, $locus"_"group[$locus];
   else
    printf gff_format, $locus, source, INTERNAL, $i, $(i+1), SCORE, $strand, remainder, $locus"_"group[$locus];
  }

  # print terminal exon
  frame=nextframe;
  remainder=($(i+1)-($i+frame)+1)%3;
  if ($strand == "+") 
    printf gff_format, $locus, source, TERMINAL, $i, $(i+1), SCORE, $strand, frame, $locus"_"group[$locus];  
  else
    printf gff_format, $locus, source, FIRST, $i, $(i+1), SCORE, $strand, remainder, $locus"_"group[$locus];  

}


