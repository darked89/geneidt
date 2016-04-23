if ($jacknifevalidate) {

jacknife_validate();
{

sub jacknife_validate {
###################################
## SET UP FOR JACKNIFE:
###################################
    #need output from processSequences4Optimization
    my $list4jkf = "";

# print STDERR "These are the sequences used for Jacknife (cross-validation) accuracy estimation\n";
    my $my_command = "grep '>' $gptrainfa | sed 's/>//g' | sort ";
    $list4jkf   = capture($my_command);

    #~ open( my $fh_LOCID, "grep '>' $gptrainfa| sed 's/>//g' | sort  |" );
    #~ while (<$fh_LOCID>) {
    #~ $list4jkf .= $_;
    #~ }
    #~ close $fh_LOCID;

    my $templist4jkf = $work_dir . $species . "_list_jacknife";
    open( my $fh_FOUT, ">", "$templist4jkf" ) or croak "Failed here";
    print $fh_FOUT "$list4jkf";
    close $fh_FOUT;

    my @list4jkf = split( "\n", $list4jkf );

    @list4jkf = sort @list4jkf;

    my $totalseqstrain = ` grep -c '^>' $gptrainfa | gawk '{print \$1}' `;
    chomp $totalseqstrain;
######
## select how many sequences to group for jacknife
######

## SET UP FOR 10 FOLD VALIDATION
    my $grpsiz = int( ( $totalseqstrain / 10 ) + 1 );

    #my $grpsiz = 1;
    print STDERR "\nThe group size for 10 fold cross validation is $grpsiz\n\n";
    my $fh_SOUT;
    print $fh_SOUT
      "\nThe group size for 10 fold cross validation is $grpsiz\n\n";

#deal with the fact that last array may have fewer elementa than the grouping value
    my $seqstofill = ( $grpsiz - @list4jkf % $grpsiz );
    if ( @list4jkf % $grpsiz ) {
        push @list4jkf, ("") x ( $grpsiz - @list4jkf % $grpsiz );
    }

    #print STDERR "$outdonorbl"."\n";
    #JACKNIFE FUNCTION

    if ( !$usebranch ) {

        @jacknifeeval = @{
            runJacknife(
                \@list4jkf,     $seqstofill,
                $grpsiz,        $outacceptortbl,
                $outdonortbl,   $outstarttbl,
                $bckgrnd,       $outcds,
                $outintron,     $bestIoWF,
                $bestIeWF,      $shortintron,
                $longintron,    $minintergenic,
                $maxintergenic, $gptraintbl,
                $gptraingff,    $startdonor,
                $enddonor,      $startacceptor,
                $endacceptor,   $startstart,
                $endstart,      0,
                0,              0,
                0,              0,
                0
            )
        };

    }
    elsif ($usebranch) {

        @jacknifeeval = @{
            runJacknife(
                \@list4jkf,          $seqstofill,    $grpsiz,
                $outacceptortbl,     $outdonortbl,   $outstarttbl,
                $bckgrnd,            $outcds,        $outintron,
                $bestIoWF,           $bestIeWF,      $shortintron,
                $longintron,         $minintergenic, $maxintergenic,
                $gptraintbl,         $gptraingff,    $startdonor,
                $enddonor,           $startacceptor, $endacceptor,
                $startstart,         $endstart,      $usebranch,
                $fullengthbranchtbl, $startbranch,   $endbranch,
                $bestAcc,            $bestMin
            )
        };

    }

    if ( !$usebranch ) {

        print STDERR
          "\n\nPerformance of new parameter file after 10x cross-validation\n\n"
          . join( "\t", @evaluationinit ), "\n";
        print $fh_SOUT
          "\n\nPerformance of new parameter file after 10x cross-validation\n\n"
          . join( "\t", @evaluationinit ), "\n";
    }
    elsif ($usebranch) {

        print STDERR
          "\n\nPerformance of new parameter file after 10x cross-validation\n\n"
          . join( "\t", @evaluationinit ), "\n";
        print $fh_SOUT
          "\n\nPerformance of new parameter file after 10x cross-validation\n\n"
          . join( "\t", @evaluationinit ), "\n";

    }

    print STDERR join( "\t", @jacknifeeval ), "\n";
    print $fh_SOUT join( "\t", @jacknifeeval ), "\n";

}    # if jacknife




if ($gff2ps) {
	gff2ps_sub(); 
}
sub gff2ps_sub {
	

#############################################################
## Obtain predictions on flanked sequences and plot them using gff2ps
#############################################################
    if ( $jacknifevalidate && $useallseqs && !$contigopt ) {

        predictPlotgff2ps(
            $paramopt,   $gptrainfa, $gptraingff,
            $gptrainlen, $temp_jkf_geneid
        );

    }
    elsif ( $contigopt && !$useallseqs ) {

        predictPlotgff2ps( $paramopt, $gpevalcontigfa, $gpevalcontiggff,
            $gpevalcontiglen, 0 );

    }
    elsif ( $contigopt && $useallseqs ) {

        predictPlotgff2ps( $paramopt, $gptraincontigfa, $gptraincontiggff,
            $gptraincontiglen, 0 );

    }
    elsif ( !$contigopt && !$useallseqs && $jacknifevalidate ) {

        predictPlotgff2ps(
            $paramopt,   $gptrainfa, $gptraingff,
            $gptrainlen, $temp_jkf_geneid
        );
    }

}

sub processBranch {

    my ( $memefile, $motifnumber, $outintron ) = @_;

#my $tempmotifbranch = $species."motif_branch";
#print STDERR "\nCHECK: begin to process branch (in function) meme: $memefile\n";
    my $motif4branch = parseMEME( $memefile, $motifnumber );

#print STDERR "\nCHECK: end of processing of branch (in function) meme: $memefile\n";
#open FOUT, ">$tempmotifbranch";
#	    print FOUT "$motif4branch";
#close FOUT;

    my $realbranches = "";
    open( my $fh_LOCID, "gawk '{print \$4}' $motif4branch  |" );
    while (<$fh_LOCID>) {
        $realbranches .= $_;
    }
    close $fh_LOCID;

    my $temprealbranches = $work_dir . $species . ".real_branches.tbl";

    open( my $fh_FOUT, ">", "$temprealbranches" ) or croak "Failed here";
    print $fh_FOUT "$realbranches";
    close $fh_FOUT;

    my $realbrancheslist = "";
    open( $fh_LOCID,
"gawk '{print \$1,\$4}' $motif4branch  | sed -e 's/^BL width.*//' -e /^\$/d |"
    );
    while (<$fh_LOCID>) {
        $realbrancheslist .= $_;
    }
    close $fh_LOCID;

    my $temprealbrancheslist = $work_dir . $species . ".real_branches_list";

    open( $fh_FOUT, ">", "$temprealbrancheslist" ) or croak "Failed here";
    print $fh_FOUT "$realbrancheslist";
    close $fh_FOUT;

    my $fullengthbranch =
      extractRealBranches( $temprealbrancheslist, $outintron );

    return $fullengthbranch;

}

sub parseMEME {

    my ( $memefi, $motifn ) = @_;
    my $motiftemp = "Motif" . $motifn;

    if ( -e "$motiftemp" ) {
        print STDERR "File $motiftemp exists already, erasing it..\n";
        unlink $motiftemp;
    }
    my $fh_MOTIF;
    my $fh_MEMEFILE;
    open( $fh_MEMEFILE, "<", "$memefi" ) || croak "Can't open $_: $!\n";
    my $in_blocks = 0;
    while (<$fh_MEMEFILE>) {
        my $line = $_;
        open( $fh_MOTIF, ">>$motiftemp" )
          || croak "Can't open $motiftemp: $!\n";
        if (/^BL\s+$fh_MOTIF $motifn/) {
            $in_blocks = 1;
        }
        if ( /\/\// && $in_blocks ) {
            $in_blocks = 0;
            close $fh_MOTIF;
        }
        elsif ($in_blocks) {
            print $fh_MOTIF "$line";
        }
    }

    #my  $fh_MEMEFILE;
    close $fh_MEMEFILE;
    close $fh_MOTIF;

    return $motiftemp;

}    #END PARSE MEME

sub extractRealBranches {
    my ( $temprealbrancheslist, $outintron ) = @_;

    my %branches;
    my %sites;

    my $fullengthbranch    = $work_dir . $species . ".real.fullbranches.tbl";
    my $fullengthbranchmod = "";

    open( my $fh_full_L_BRANCH, ">", "$fullengthbranch" )
      or croak "Failed here";
    open( my $fh_BRANCHES, "<", "$temprealbrancheslist" )
      || croak "Can't open branches\n";

    while (<$fh_BRANCHES>) {
        my ( $id, $sequence ) = split;
        $branches{$id} = $sequence;
    }

    close($fh_BRANCHES);

    #??XX close

    open( my $fh_INTRONS, "<", "$outintron" ) || croak "Can't open introns\n";

    my $bpindex = $bp - 1;
    while (<$fh_INTRONS>) {
        my ( $id, $sequence ) = split;
        $sequence = sprintf "%s%s%s", 'N' x 30, $sequence, 'N' x 30;
        foreach my $key (%branches) {
            if (   ( $id =~ /$key/ )
                && ( my $pos = index( $sequence, $branches{$key} ) ) )
            {
                $sites{$key} = substr( $sequence, $pos + $bpindex - 32, 62 )
                  ;   ### ORIGINALLY: $sites{$key}=substr($sequence,$pos-21,62);
                print $fh_full_L_BRANCH "$id\t$sites{$key}\n";
            }
        }
    }
    close($fh_BRANCHES);
    close($fh_full_L_BRANCH);

    my $fullengthbranches =
      $work_dir . $species . ".real.full_length_branches.tbl";
    open( my $fh_FOUT, ">", "$fullengthbranches" ) or croak "Failed here";
    print $fh_FOUT "$fullengthbranchmod";
    close $fh_FOUT;

    unlink $fullengthbranch;

    return $fullengthbranches;

}

## JACKNIFE SUB
sub runJacknife {

    my (
        $list4jkf,           $seqstofill,    $grpsiz,
        $outacceptortbl,     $outdonortbl,   $outstarttbl,
        $bckgrnd,            $outcds,        $outintron,
        $bestIoWF,           $bestIeWF,      $shortintron,
        $longintron,         $minintergenic, $maxintergenic,
        $gptraintbl,         $gptraingff,    $startdonor,
        $enddonor,           $startacceptor, $endacceptor,
        $startstart,         $endstart,      $branchsw,
        $fullengthbranchtbl, $startbranch,   $endbranch,
        $bestAcc,            $bestMin
    ) = @_;

    #my $grepBranches = "";
    my $branchexcept = "";
    my $count        = 1;
    my $geneidout    = "";

    #my $cycles = 402;
    my $cycles = int( scalar(@$list4jkf) / $grpsiz );    #610/10
    my $status = "";

    #  my $tempgroup4jckf = "";
    my $actualannots = scalar(@$list4jkf) - $seqstofill;

    print STDERR
"For the 10-fold cross-validation method we will remove groups of $grpsiz sequences from the "
      . $actualannots
      . " training set annotations, train on the remaining annotations and run geneid evaluate the performance on the left-out sequences (performed $cycles times)\n\n";
    print $fh_SOUT
"For the 10-fold cross-validation method we will remove groups of $grpsiz sequences from the "
      . $actualannots
      . " training set annotations, train on the remaining annotations and run geneid evaluate the performance on the left-out sequences (performed $cycles times)\n\n";

    ## for loop to go through seqs
    for ( my $i = 0 ; $i < scalar(@$list4jkf) ; $i += $grpsiz ) {

        print STDERR
"Extracting the set $count of sequences for 10x cross validation training and evaluation ($count out of $cycles)\n";

        my @seqsjk = @$list4jkf[ $i .. $i + $grpsiz - 1 ];
        my $listseqsjacknife = join( "\n", @seqsjk );

        #        print STDERR "group1:".$listseqsjacknife."\n\n";
        # my $listseqsjacknife4eval =  "$listseqsjacknife"."#$";

        my $group4jckf = $work_dir . $species . "_seqs_leftout";
        open( my $fh_FOUT, ">", "$group4jckf" ) or croak "Failed here";
        print $fh_FOUT "$listseqsjacknife";
        close $fh_FOUT;

##grepSeqsToEval

        my $seqsaside1 = "";
        open( my $fh_LOCID,
            "sed -e '/^\$/d' -e 's/\\(.*\\)/\\1\$/g' $group4jckf |" );
        while (<$fh_LOCID>) {
            $seqsaside1 .= $_;
        }
        close $fh_LOCID;

        my $tempgroup4jckf1 =
          $work_dir . $species . "_seqs_leftout_crossvalidation_1";
        open( $fh_FOUT, ">", "$tempgroup4jckf1" );
        print $fh_FOUT "$seqsaside1";
        close $fh_FOUT;

        #		print STDERR " $tempgroup4jckf1\n";
        #Acceptors,Donors,Starts, Noncoding,branches

        my $seqsaside2 = "";
        open( $fh_LOCID,
            "sed -e '/^\$/d' -e 's/\\(.*\\)/\\1\./g' $group4jckf |" );
        while (<$fh_LOCID>) {
            $seqsaside2 .= $_;
        }
        close $fh_LOCID;

        my $tempgroup4jckf2 =
          $work_dir . $species . "_seqs_leftout_crossvalidation_2";
        open( $fh_FOUT, ">", "$tempgroup4jckf2" );
        print $fh_FOUT "$seqsaside2";
        close $fh_FOUT;

        #	         print STDERR " $tempgroup4jckf2\n";

        #Coding

        my $seqsaside3 = "";
        open( $fh_LOCID,
            "sed -e '/^\$/d' -e 's/\\(.*\\)/\\1\_/g' $group4jckf |" );
        while (<$fh_LOCID>) {
            $seqsaside3 .= $_;
        }
        close $fh_LOCID;

        my $tempgroup4jckf3 =
          $work_dir . $species . "_seqs_leftout_crossvalidation_3";
        open( $fh_FOUT, ">", "$tempgroup4jckf3" );
        print $fh_FOUT "$seqsaside3";
        close $fh_FOUT;

        #	                 print STDERR " $tempgroup4jckf3\n";

        #CREATE BLANK TEMP PARAMETER FILE FOR JACKNIFE
        my $paramtemp = Geneid::Param->new($species);

        #set isochores to 1
        $paramtemp->numIsocores(1);
        $paramtemp->isocores( [ Geneid::Isocore->new() ] );

        #	print SDTERR $outacceptortbl."\n";
## PRE EXTRACT SEQUENCES FROM TRAINING SET
        my $my_command;
        $my_command =
"egrep -vf $tempgroup4jckf2 $outacceptortbl > $tmp_dir/tmp_Acceptorsminus";
        run($my_command);
        $my_command =
"egrep -vf $tempgroup4jckf2 $outdonortbl    > $tmp_dir/tmp_Donorsminus";
        run($my_command);
        $my_command =
"egrep -vf $tempgroup4jckf2 $outstarttbl    > $tmp_dir/tmp_Startsminus";
        run($my_command);
        $my_command =
"egrep -vf $tempgroup4jckf3 $outcds         > $tmp_dir/tmp_Codingminus";
        run($my_command);
        $my_command =
"egrep -vf $tempgroup4jckf2 $outintron      > $tmp_dir/tmp_NonCodingminus";
        run($my_command);

        $my_command =
"gawk '{print \$2,\$1}' $gptraintbl | egrep -f $tempgroup4jckf1 - | gawk '{print \$2,\$1}' - > $tmp_dir/tmp_SeqsToEval";
        run($my_command);

        #~ #	print STDERR "$gptraintbl $tempgroup4jckf1";

#~ my $grepSeqsToEval
#~ = "gawk '{print \$2,\$1}' $gptraintbl | egrep -f $tempgroup4jckf1 - | gawk '{print \$2,\$1}' - > $tmp_dir/tmp_SeqsToEval";

        if ($branchsw) {
            $my_command =
"egrep -vf $tempgroup4jckf2 $fullengthbranchtbl  > $tmp_dir/tmp_Branchesminus";
            run($my_command);

#$grepBranches
#    = "egrep -vf $tempgroup4jckf2 $fullengthbranchtbl  > $tmp_dir/tmp_Branchesminus";

        }

## store taken out sequences in the following variables:

        my $acceptorexcept  = "$tmp_dir/tmp_Acceptorsminus";
        my $donorexcept     = "$tmp_dir/tmp_Donorsminus";
        my $startexcept     = "$tmp_dir/tmp_Startsminus";
        my $codingexcept    = "$tmp_dir/tmp_Codingminus";
        my $noncodingexcept = "$tmp_dir/tmp_NonCodingminus";
        my $seqstoevaltbl   = "$tmp_dir/tmp_SeqsToEval";
        if ($branchsw) {
            $branchexcept = "$tmp_dir/tmp_Branchesminus";
        }

## CONVERT TO FASTA (EVALUATE SET)
        my $tempgp_fa_eval_jkf =
          $work_dir . $species . ".gp.eval.crossvalidation.fa";
        TblToFasta( $seqstoevaltbl, $tempgp_fa_eval_jkf );

## EXTRACT SEQUENCES FROM JACKNIFE TRAINING SET

#########
## get donor site statistics
#########

        my $order = "0";
        $numbersites = num_of_lines_in_file($donorexcept);

        #my $numbersites = ` wc -l $donorexcept | gawk '{print \$1}'`;
        #chomp $numbersites;
        #$numbersites = int($numbersites);
        my $donoffset = $bases_offset
          ;    #position before intron (last of exon (31) -1 for offset)

        if ( $numbersites > $train_sites_cutoff_alt ) {

            $order = "1";

        }
        elsif ( $numbersites <= $train_sites_cutoff_alt ) {

            $order = "0";
        }

        print STDERR
"$startdonor/$enddonor/$prof_len_don There are $numbersites donor sites, enough for a matrix of order $order (Jacknife) set $count \n";

        my ( $donormatrix, $prof_len_don, $fxddonoffset, $s, $e ) =
          getKmatrix( $donorexcept, $bckgrnd, $order, $donoffset, 1, 0, 0, 0,
            $startdonor, $enddonor, 1 );
        if (
            !defined @{ $paramtemp->isocores }[0]->set_profile(
                'Donor_profile', $prof_len_don, $fxddonoffset, $pwm_cutoff,
                $order, 0, 1, 0, 0, 0, 0, $donormatrix
            )
          )
        {
            croak "error in setting profile\n";
        }

########($true_seqs,$false_seqs,$order,$offset,$donor,$accept,$star,$branch,$start,$end,$jacknife)

## get acceptor site statistics

        $order       = "0";
        $numbersites = num_of_lines_in_file($acceptorexcept);

        #$numbersites = ` wc -l $acceptorexcept | gawk '{print \$1}'`;
        #chomp $numbersites;
        #$numbersites = int($numbersites);
        my $accoffset = $bases_offset
          ;    #position after intron (first of exon (31) -1 for offset)

        if ( $numbersites > $train_sites_cutoff_alt ) {

            $order = "1";

        }
        elsif ( $numbersites <= $train_sites_cutoff_alt ) {

            $order = "0";
        }

        print STDERR
"$startacceptor/$endacceptor/$prof_len_acc:There are $numbersites acceptor sites, enough for a matrix of order $order (Jacknife) set $count \n";

        my ( $acceptormatrix, $prof_len_acc, $fxdaccoffset, $st, $en ) =
          getKmatrix( $acceptorexcept, $bckgrnd, $order, $accoffset, 0, 1, 0,
            0, $startacceptor, $endacceptor, 1 );
        if (
            !defined @{ $paramtemp->isocores }[0]->set_profile(
                'Acceptor_profile', $prof_len_acc, $fxdaccoffset, $pwm_cutoff,
                $order, 0, 1, 0, 0, 0, 0, $acceptormatrix
            )
          )
        {
            croak "error in setting profile\n";
        }

#########my ($true_seqs,$false_seqs,$order,$offset,$donor,$accept,$star,$branch,$start,$end,$jacknife)

## get start site statistics

        $order       = "0";
        $numbersites = num_of_lines_in_file($startexcept);

        #$numbersites = ` wc -l $startexcept | gawk '{print \$1}'`;
        #chomp $numbersites;
        #$numbersites = int($numbersites);
        my $staoffset = $bases_offset
          ;    #before first position of the exon (31)minus 1 for offset)

        if ( $numbersites > $train_sites_markov_cutoff ) {

            $order = "2";

        }
        elsif ( $numbersites <= $train_sites_markov_cutoff ) {

            $order = "0";
        }

        print STDERR
"There are $numbersites start sites, enough for a matrix of order $order (Jacknife) set $count \n";

        my ( $startmatrix, $prof_len_sta, $fxdstaoffset, $sta, $enj ) =
          getKmatrix( $startexcept, $bckgrnd, $order, $staoffset, 0, 0, 1, 0,
            $startstart, $endstart, 1 );

####write to parameter file
        if (
            !defined @{ $paramtemp->isocores }[0]->set_profile(
                'Start_profile', $prof_len_sta, $fxdstaoffset, $pwm_cutoff,
                $order, 0, 1, 0, 0, 0, 0, $startmatrix
            )
          )
        {
            croak "error in setting profile\n";
        }
#############################my ($true_seqs,$false_seqs,$order,$offset,$donor,$accept,$star,$branch,$start,$end,$jacknife) = @_;

        if ($branchsw) {

## get branch site statistics

            $order       = "0";
            $numbersites = num_of_lines_in_file($branchexcept);

            #$numbersites = ` wc -l $branchexcept | gawk '{print \$1}'`;
            #chomp $numbersites;
            #$numbersites = int($numbersites);
            my $braoffset =
              "32";   #before first position of the exon (35)minus 1 for offset)

            print STDERR
"There are $numbersites branch sites, enough for a matrix of order $order (Jacknife) set $count \n";

            my (
                $branchmatrix, $prof_len_bra, $fxdbraoffset,
                $startbranch,  $endbranch
              )
              = getKmatrix( $fullengthbranchtbl, $bckgrnd, $order,
                $braoffset, 0, 0, 0, 1, $startbranch, $endbranch, 1 );

## write to parameter file
            if (
                !defined @{ $paramtemp->isocores }[0]->set_profile(
                    'Branch_point_profile', $prof_len_bra, $fxdbraoffset, "-50",
                    $order, 0, 1, $bestAcc, $bestMin, 0, 0, $branchmatrix
                )
              )
            {
                croak "error in setting profile\n";
            }
#############################

        }

## DERIVE INITIAL/TRANSITION MARKOV MODEL JACKNIFE

       #print STDERR "\nDeriving markov model (Jacknife method): set $count \n";

        my (
            $markovini,      $markovtrans, $totalcoding,
            $totalnoncoding, $markovmodel
        ) = @{ deriveCodingPotential( $codingexcept, $noncodingexcept ) };

        #add markov matrices to the parameter file
        if ( !defined @{ $paramtemp->isocores }[0]->Markov_order($markovmodel) )
        {
            croak "error in setting Markov_order\n";
        }
        if ( !defined @{ $paramtemp->isocores }[0]
            ->Markov_Initial_probability_matrix($markovini) )
        {
            croak "error in setting Markov_Initial_probability_matrix\n";
        }
        if ( !defined @{ $paramtemp->isocores }[0]
            ->Markov_Transition_probability_matrix($markovtrans) )
        {
            croak "error in setting Markov_Transition_probability_matrix\n";
        }

## WRITE PRELIMINARY PARAMETER FILE FOR JACKNIFE
        #$paramtemp->writeParam("$species.geneid.param.tmp");

        #my $newparamtemp = "$species.geneid.param.tmp";

        #$paramtemp = Geneid::Param->new();
        #	     $param->readParam("$newparamtemp");

        for ( my $i = 0 ; $i < $paramtemp->numIsocores ; $i++ ) {
            if (
                !defined @{ $paramtemp->isocores }[$i]->Exon_weights(
                    [ $bestIeWF, $bestIeWF, $bestIeWF, $bestIeWF ]
                )
              )
            {
                croak "error in setting exon weights\n";
            }
            if ( !defined @{ $paramtemp->isocores }[$i]
                ->Exon_factor( [ $bestIoWF, $bestIoWF, $bestIoWF, $bestIoWF ] )
              )
            {
#   if (!defined @{$paramtemp->isocores}[$i]->Exon_factor([0.4,$bestIoWF,$bestIoWF,0.4])) {
                croak "error in setting exon weights\n";
            }
            if (
                !defined @{ $paramtemp->isocores }[$i]->Site_factor(
                    [
                        1 - $bestIoWF,
                        1 - $bestIoWF,
                        1 - $bestIoWF,
                        1 - $bestIoWF
                    ]
                )
              )
            {
#  if (!defined @{$paramtemp->isocores}[$i]->Site_factor([0.55,1-$bestIoWF,1-$bestIoWF,0.55])) {
                croak "error in setting exon weights\n";
            }
        }    # for number of isochores

        #print STDERR "bestIoWF jacknife: $bestIoWF\n";
        #Open gene model object
        $paramtemp->geneModel( Geneid::GeneModel->new() );
        $paramtemp->geneModel->useDefault;
        $paramtemp->geneModel->intronRange( $shortintron, $longintron );

        #$paramtemp->geneModel->intergenicRange(200,'Infinity');
        $paramtemp->geneModel->intergenicRange( $minintergenic,
            $maxintergenic );
## WRITE PRELIMINARY PARAMETER FILE FOR JACKNIFE
        $paramtemp->writeParam("$species.geneid.jacknife.param.temp.$count");
        my $newparamtemp = "$species.geneid.jacknife.param.temp.$count";

## RUN GENEID ON THE LEFT OUT SEQUENCES:
        print STDERR "\nRunning geneid on set number $count\n\n\n";

        open $fh_LOCID,
"./bin/geneid -GP $newparamtemp $tempgp_fa_eval_jkf | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' | ";
        while (<$fh_LOCID>) {
            $geneidout .= $_;
        }
        close $fh_LOCID;

        $geneidout .= '#$' . "\n";

        $count++;

    }    #end of for loop				###FOR LOOP: GO THROUGH GROUP OF SEQS

## EVALUATE GENEID AGAINST ANNOTATIONS

    $temp_jkf_geneid = $work_dir . $species . ".geneid_crossvalidation";
    open( my $fh_FOUT, ">", "$temp_jkf_geneid" );
    print $fh_FOUT "$geneidout";
    close $fh_FOUT;

    if ($branchsw) {

        my @evaluation_jacknife = split " ",
` ./bin/evaluation -sta $temp_jkf_geneid $gptraingff | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\%6.2f \%6.2f \%6.2f \%6.2f\\n\", $bestIoWF,$bestIeWF,$bestAcc,$bestMin,\$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}' `;

        return \@evaluation_jacknife;

    }
    elsif ( !$branchsw ) {

        my @evaluation_jacknife = split " ",
` ./bin/evaluation -sta $temp_jkf_geneid $gptraingff | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $bestIoWF,$bestIeWF,\$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}' `;

        return \@evaluation_jacknife;

    }
    return 1;

    #unlink $temp_jkf_geneid;

}    #SUB JACKNIFE

