#!/usr/bin/gawk -f  
# rgs, imim.sb, aug 98
# gff2cds convers a gff file to a cds file
# 
BEGIN{
  # defaults

  USAGE="USAGE: gff2cds source=\"source\"";

  # gff data structure
  seqname  = 1;
  source_  = 2;
  feature  = 3;
  start    = 4;
  end      = 5;
  score    = 6;
  strand   = 7;
  frame    = 8;
  group    = 9;

}
{
  if (!source) {
    print USAGE;
    error=1; exit;
  }
}    

$source_==source{

  INDEX0=$seqname SUBSEP $group;

  # put cds in the correct position --sorted by start position
  for (j=1;j<=ncds[INDEX0] && $start > cds[INDEX0,j,1]; j++);
  
  for (k=ncds[INDEX0];k>=j;k--) {
    cds[INDEX0,k+1,1]=cds[INDEX0,k,1];
    cds[INDEX0,k+1,2]=cds[INDEX0,k,2];
  }
  cds[INDEX0,j,1]=$start;
  cds[INDEX0,j,2]=$end;
  ncds[INDEX0]++;
  strand_[INDEX0]=$strand;
  seqname_[$group]=$seqname;
}
END{
  if (error)
    exit;
  for (g in ncds) {
    split (g,gx,SUBSEP);
    printf "%s %s", gx[1], strand_[g];
    for (i=1;i<=ncds[g];i++)
      printf  " %d %d", cds[g,i,1], cds[g,i,2];
    printf "\n";
  }
}
