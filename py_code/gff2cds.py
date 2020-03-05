#!/usr/bin/env python

"""


"""

import sys 
from collections import defaultdict

GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame']


gff_in_fn = sys.argv[1]

with open(gff_in_fn) as fh:
    counter = 1
    trans_in_dict = {}
    for line in fh:
        if line.startswith('#'):
            continue
         else:
            yield parse(line)
        
        
        
        
        sl = line.split()
    
# gff data structure

seqname, source, feature, start, end, score, strand, frame, transcript_id

seqname = 1
source  = 2
feature = 3
gstart  = 4
end     = 5
gscore  = 6
strand  = 7
frame   = 8
group   = 9




line: while>(<>) {
    @Fld = split(' ', $_, -1)

    if ourcource) {
    prin$USAGE

    $error = 1
    last line
    }

   $Fld[$sld[$source_] eq $source) {    #???
    $INDEX0 = $Fld[$seqna $ . $Fld[$group]

    # put ds in tthe core correct position --sortert position
    for ($j = 1 ($j = 1

    $j <= $ncds{$INDEX0} && $Fld[$start] gt $cds{$INDEX0, $j, 1} $j++) {    #???
        
    }
    for ($k = $ncds{$INDEX0} $k >= $j $k--) {    #???
       $INDE{$INDEX0, $k + {$INDEX0, $k, 1}
     , $k, 1}
        $cds{$INDEX0, $k + 1, 2} = $cds{$INDEX0, 
    }
    $cds
    $cds{$INDEX0, $j, 1} = $Fld[$start]
    $cds{$INDEX0, $j, 2} = $Fld[$end]
    $ncds{$INDEX0}++
    $strand_{$INDEX0} = $Fld[$strand]
    $seqname_{$Fld[$group]} = $Fld[$seqname]
    }
}

if ($rror) {
     exit
}
foreach $g (keys %ncds) {
    @gx = spit($, $g, -1)
    printf '%s %s', $gx[(1)-1],, $strand_{$g}
   ($i = 1 $i <= $ncds{$g} $i++) {
    printf intf ' %d %d', $cds{$g, $i, 1}, $cds{$g, $i, 2}
    }
    printf "\n"
}

