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


sub start_reduced {
#############################################################
## reduced/short version training starting with PWMs PLUS BACKGROUND#
#############################################################
    ###reduced/short training starting with PWMs PLUS BACKGROUND
    $reducedtraining = 1;
    print STDERR
"\n\nYou chose to continue the training of $species by assuming the cds, intronic sequences and splice sites have already been extracted... \n\n";
## XX BUG: we did not store the used variables + subsequent
## XX BUG: eval $_ is  against the good PERL practice

#open( my $fh_species_VARS, "<", "${species}.variables" )
#|| die
#"You need to have run the training program once previously to execute the reduced version of the geneid training program \n";
#while (<$fh_species_VARS>) {
#eval $_;
#}
#die "can't restore variables from ${species}.variables: $@" if $@;
#close $fh_species_VARS;

## CREATE A STATS FILE
    my @timeData = localtime(time);

    #STATS DIR CREATED FIRST TIME PIPELINE IS RUN FOR A GIVEN SPECIES
    my $statsout = $stats_dir . join( '_', @timeData ) . "_$sout";
## OPEN STATISTICS OUTPUT AT THIS TIME...EVERY TIME PIPELINE IS RUN
    open( my $fh_SOUT, ">", "$statsout" ) or croak "Failed here";

    if ( !$useallseqs ) {
        print STDERR
"\nThe reduced training process will use 80% of the gene-model sequences ($totalseqs4training)/20% will used for posterior evaluation of the newly developed parameter file ($gffseqseval)\n";
    }
    else {
        print STDERR
"The reduced training process will use ALL of the gene-model sequences ($total_seqs)\n";
    }
    if ($jacknifevalidate) {
        print STDERR
"\nThe reduced training process will include a 10x cross validation of the accuracy of the new $species parameter file\n";
    }

    print STDERR
"\nA subset of $totalseqs4training sequences (randomly chosen from the $total_seqs gene models) was used for training\n";
    print $fh_SOUT "GENE MODEL STATISTICS FOR $species\n\n";
    return 1;
}    #end reduced


sub refactor_branch_sub {

        #print STDERR "\nCHECK: begin to process branch\n";
        $fullengthbranchtbl =
          processBranch( $memefile, $motifnumber, $outintron );

#########
## get branch site statistics
#########

        $order = "0";

        $numbersites = num_of_lines_in_file($fullengthbranchtbl);

        my $braoffset = "32"
          ; #before the A (branch) (33)minus 1 for offset) ####CHANGE TO 26 for it to work with meme.txt of V. dahliae introns...

        print STDERR
"\nThere are $numbersites branch sites, enough for a matrix of order $order, offset: $braoffset \n";

        (
            $branchmatrix, $prof_len_bra, $fxdbraoffset,
            $startbranch,  $endbranch
          )
          = getKmatrix( $fullengthbranchtbl, $bckgrnd, $order, $braoffset, 0,
            0, 0, 1, 0, 0, 0 );

## write to parameter file
        if (
            !defined @{ $param->isocores }[0]->set_profile(
                'Branch_point_profile', $prof_len_bra, $fxdbraoffset, -50,
                $order, 0, 1, 40, 10, 0, 0, $branchmatrix
            )
          )
        {
            croak "error in setting profile\n";
        }

        my $brasub = "";

#print STDERR "gawk '{print  substr(\$2,($startstart-3),($prof_len_don+6))}' $outstarttbl\n";
        $my_command =
"gawk '{print  substr(\$2,($startbranch-3),($prof_len_bra+6))}' $fullengthbranchtbl ";
        $brasub = capture($my_command);


        $branchsubprofile = $work_dir . $species . ".bra.sub.profile";
        open( my $fh_FOUT, ">", "$branchsubprofile" ) or croak "Failed here";
        print $fh_FOUT "$brasub";
        close $fh_FOUT;

 #print STDERR "$path/pictogram $startsubprofile $statsdir/Start -bits -land\n";
## BUG?
#    run("./bin/pictogram $branchsubprofile ./statistics_${species}/Branch -bits -land");
        $my_command =
"./bin/pictogram $branchsubprofile $plots_dir/branch_profile.pictogram -bits -land";
        print "\n$my_command\n";
        run($my_command);

        #unlink $branchsubprofile;
        return 1;

    } #  refactor use branch


sub predictPlotgff2ps {
    my ( $paramopt, $gpfa, $gpgff, $gplen, $tempjkf_geneid ) = @_;

    print STDERR " \nRunning geneid on fastas \n \n ";

    my $geneidall   = "";
    my $gff2psplots = "$species . gff2ps . prediction . plots . ps ";
    ## TOFIX
    my $my_command =
" . /bin/geneid -GP $paramopt $gpfa | gawk
'NR>5 {OFS=\"\\t\";if (\$3==\"Gene\") print \"\#\$\"; \$2=\"geneid_$species\"; if (substr(\$1,1,1)!=\"\#\") print }'
                              | egrep -wv 'exon' ";
    $geneidall = capture($my_command);

    my $tempgeneidgffpreds = $work_dir . $species . "
                              . geneid
                              . predictions
                              . gff ";
    open( my $fh_FOUT, " > ", "$tempgeneidgffpreds " ) or croak " Failed here ";
    print $fh_FOUT "$geneidall ";
    close $fh_FOUT;

    unlink $gff2psplots;

    #print STDERR " Plotting of predictions using gff2ps \n ";

    open( my $fh_LOCI_gplen, " < ", "$gplen " ) or croak " Failed here ";
    print STDERR " \nCreating gff2ps plots for $species \n \n ";
    while (<$fh_LOCI_gplen>) {
        my ( $gene_id, $genelength ) = split;
        print STDERR " \n $gplen $gpgff $tempjkf_geneid $gene_id $genelength \n
                              ";
        run(" egrep -w '$gene_id' $gpgff > $tmp_dir$gene_id . gff ");
        run(" egrep -w '$gene_id' $tempgeneidgffpreds >> $tmp_dir$gene_id
                              . gff ");
        if ($jacknifevalidate) {
            run(" egrep -w '$gene_id' $tempjkf_geneid >> $tmp_dir$gene_id
                              . gff ");
        }
        if ( !$contigopt ) {
            run(
"
                              . /bin/gff2ps -v -p --$tmp_dir$gene_id
                              . gff > $plots_dir / $species
                              . ${gene_id}
                              . ps "
            );
            print STDERR "    #";
                          } elsif ($contigopt) {
                            my $nucleotidesperline = 10000;
                            run(
" ./bin/gff2ps -v -p -N $nucleotidesperline -C $work_dir.gff2psrcNEW -- $tmp_dir$gene_id.gff > $plots_dir/$species.gv"
                            );
                            print STDERR "#";
                        }
                    }
                    close $fh_LOCI_gplen;
                    return 1;
}


## BUG unused. we do not use interactive or branch in simplified script 
#~ elsif ($usebranch) {    #use separate branch profile
    #~ if ($interactive) {
        #~ my $respo = "";
        #~ do {
            #~ print STDERR
#~ "Use automatically selected range values for the optimization of geneid eWF (exon weight)/oWF (exon/oligo factor)/Minimum Branch Distance from Acceptor and Context length in which Branch sites should be scored (AccCtx) internal parameters?\n\n(eWF: $IeWF to $FeWF; step $deWF\noWF: $IoWF to $FoWF; step $doWF\nMinBranchDistance: $iMin to $fMin; step $dMin\nAcceptorBranchCtx: $iAccCtx to $fAccCtx; step $dAccCtx)\n\nDo you prefer to change these values? (we do not recommed you change the Min Branch Distance and AccBranch context)";
            #~ $respo = readline(STDIN);
        #~ } while ( $respo !~ /^(yes|y)|(n|no)$/i );

        #~ if ( $respo =~ /^(yes|y)/i ) {

            #~ my $sline = "";
            #~ my $eline = "";
            #~ my $dline = "";
            #~ do {
                #~ print STDERR "\nType new initial eWF (IeWF): ";
                #~ $sline = readline(STDIN);
            #~ } while ( $sline !~ /(-*[0-9]*\.*[0-9]+)/ );
            #~ $IeWF = $1;

            #~ do {
                #~ print STDERR "\nType new final eWF (FeWF): ";
                #~ $eline = readline(STDIN);
              #~ } while ( $eline !~ /(-*[0-9]*\.*[0-9]+)/
                #~ || $eline <= $sline );
            #~ $FeWF = $1;

            #~ do {
                #~ print STDERR "\nType step (delta) eWF (deWF)): ";
                #~ $dline = readline(STDIN);
            #~ } while ( $dline !~ /(-*[0-9]*\.*[0-9]+)/ );
            #~ $deWF = $1;

            #~ do {
                #~ print STDERR "\nType new initial oWF (IoWF): ";
                #~ $sline = readline(STDIN);
            #~ } while ( $sline !~ /(-*[0-9]*\.*[0-9]+)/ );
            #~ $IoWF = $1;

            #~ do {
                #~ print STDERR "\nType new final oWF (FoWF): ";
                #~ $eline = readline(STDIN);
              #~ } while ( $eline !~ /(-*[0-9]*\.*[0-9]+)/
                #~ || $eline <= $sline );
            #~ $FoWF = $1;

            #~ do {
                #~ print STDERR "\nType step (delta) oWF (doWF): ";
                #~ $dline = readline(STDIN);
            #~ } while ( $dline !~ /(-*[0-9]*\.*[0-9]+)/ );
            #~ $doWF = $1;

            #~ do {
                #~ print STDERR "\nType new initial Min Branch Distance (iMin): ";
                #~ $sline = readline(STDIN);
            #~ } while ( $sline !~ /(-*[0-9]*\.*[0-9]+)/ );
            #~ $iMin = $1;

            #~ do {
                #~ print STDERR "\nType new final Min Branch Distance (fMin): ";
                #~ $eline = readline(STDIN);
              #~ } while ( $eline !~ /(-*[0-9]*\.*[0-9]+)/
                #~ || $eline <= $sline );
            #~ $fMin = $1;

            #~ do {
                #~ print STDERR
                  #~ "\nType step (delta) Min Branch Distance (dMin)): ";
                #~ $dline = readline(STDIN);
            #~ } while ( $dline !~ /(-*[0-9]*\.*[0-9]+)/ );
            #~ $dMin = $1;

            #~ do {
                #~ print STDERR
                  #~ "\nType new initial Acceptor/Branch Context (iAccCtx): ";
                #~ $sline = readline(STDIN);
            #~ } while ( $sline !~ /(-*[0-9]*\.*[0-9]+)/ );
            #~ $iAccCtx = $1;

            #~ do {
                #~ print STDERR
                  #~ "\nType new final Acceptor/Branch Context (fAccCtx): ";
                #~ $eline = readline(STDIN);
              #~ } while ( $eline !~ /(-*[0-9]*\.*[0-9]+)/
                #~ || $eline <= $sline );
            #~ $fAccCtx = $1;

            #~ do {
                #~ print STDERR
                  #~ "\nType step (delta) Acceptor/Branch Context (dAccCtx): ";
                #~ $dline = readline(STDIN);
            #~ } while ( $dline !~ /(-*[0-9]*\.*[0-9]+)/ );
            #~ $dAccCtx = $1;

        #~ }
    #~ }
#~ ## OPTIMIZATION FUNCTIONS


sub branch_start {
    ###IF THERE IS A MEME-DISCOVERED BRANCH POINT PROFILE
    $usebranch = 1;
    do {
        print STDERR
"\n\nYou chose to use meme branch profile information in training $species. Please, indicate the name of the meme output file (make sure the file is in the right path) followed by the number of the motif the branch profile and the position(not index) of the branchpoint A (i.e. enter 'meme.txt 2 7' if the meme output file is called meme.txt and the second motif is the one corresponding to the branch profile and the A is at the seventh position)\n";
        $memeanswer = <ARGV>;
        chomp $memeanswer;
    } while ( $memeanswer !~ /^(.+)\s+(\d+)\s+(\d+)$/i );

    $memefile    = $1;
    $motifnumber = $2;
    $bp          = $3;

    if ( -e "$memefile" ) {
        print "File exists and is named: ($memefile) \n\n";
    }
    else {
        print "File does not exist";

        print STDERR $usage and exit;

    }
    return 1;
}



## not complete function for branch opty
sub OptimizeParameter {

    my (
        $gpfa,         $gpgff,        $newparam,     $branchswitch,
        $prof_len_bra, $fxdbraoffset, $branchmatrix, $IeWF,
        $deWF,         $FeWF,         $IoWF,         $doWF,
        $FoWF,         $iMin,         $dMin,         $fMin,
        $iAccCtx,      $dAccCtx,      $fAccCtx
    ) = @_;
    my @evaluation_total = ();
    my $IeWFini          = $IeWF;
    my $IoWFini          = $IoWF;
    my $iMinini          = $iMin;
    my $iAccCtxini       = $iAccCtx;

    my $fh_SOUT;

    open( $fh_SOUT, ">", "$work_dir/$species.OptimizeParameter.log" )
      or croak "Failed here";

    #~ if ( !$branchswitch ) {
        
        print STDERR
          "\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\n\n";
        print $fh_SOUT
          "\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\n\n";
    #~ }
    #~ else {
        #~ print STDERR
#~ "\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\nMinBranch range : $iMin to $fMin\nAccCtx range : $iAccCtx to $fAccCtx\n\n";
        #~ print $fh_SOUT
#~ "\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\nMinBranch range : $iMin to $fMin\nAccCtx range : $iAccCtx to $fAccCtx\n\n";

    #~ }

    for ( $IeWF = $IeWFini ; $IeWF <= $FeWF ; $IeWF += $deWF ) {
        print STDERR "eWF: $IeWF\noWF: ";

        for ( $IoWF = $IoWFini ; $IoWF <= $FoWF ; $IoWF += $doWF ) {
			print STDERR "$IoWF  ";
			
            #~ if ( !$branchswitch ) {
                #~ print STDERR "$IoWF  ";
            #~ }
            #~ if ($branchswitch) {
                #~ print STDERR "$IoWF\nMinDist:  ";
                #~ for ( $iMin = $iMinini ; $iMin <= $fMin ; $iMin += $dMin ) {
                    #~ print STDERR "$iMin\nAccCtx:  ";

                    #~ for (
                        #~ $iAccCtx = $iAccCtxini ;
                        #~ $iAccCtx <= $fAccCtx ;
                        #~ $iAccCtx += $dAccCtx
                      #~ )
                    #~ {
                        #~ print STDERR "$iAccCtx  ";
                        #~ my $param = Geneid::Param->new();
                        #~ $param->readParam("$newparam");

                        #~ for ( my $i = 0 ; $i < $param->numIsocores ; $i++ ) {
                            #~ if (
                                #~ !defined @{ $param->isocores }[$i]
                                #~ ->Exon_weights(
                                    #~ [ $IeWF, $IeWF, $IeWF, $IeWF ]
                                #~ )
                              #~ )
                            #~ {
                                #~ croak "error in setting exon weights\n";
                            #~ }
                            #~ if ( !defined @{ $param->isocores }[$i]
                                #~ ->Exon_factor( [ $IoWF, $IoWF, $IoWF, $IoWF ] )
                              #~ )
                            #~ {
#~ # if (!defined @{$param->isocores}[$i]->Exon_factor([0.33,$bestIoWF,$bestIoWF,0.33])) {
                                #~ croak "error in setting exon weights\n";
                            #~ }
                            #~ if (
                                #~ !defined @{ $param->isocores }[$i]
                                #~ ->Site_factor(
                                    #~ [
                                        #~ 1 - $IoWF, 1 - $IoWF,
                                        #~ 1 - $IoWF, 1 - $IoWF
                                    #~ ]
                                #~ )
                              #~ )
                            #~ {
#~ # if (!defined @{$param->isocores}[$i]->Site_factor([0.45,1-$bestIoWF,1-$bestIoWF,0.45])) {
                                #~ croak "error in setting exon weights\n";
                            #~ }
                            #~ if (
                                #~ !defined @{ $param->isocores }[$i]
                                #~ ->set_profile(
                                    #~ 'Branch_point_profile', $prof_len_bra,
                                    #~ $fxdbraoffset,          -50,
                                    #~ 0,                      0,
                                    #~ 1,                      $iAccCtx,
                                    #~ $iMin,                  0,
                                    #~ 0,                      $branchmatrix
                                #~ )
                              #~ )
                            #~ {
                                #~ croak "error in setting profile\n";
                            #~ }
                        #~ }
                        #~ my $temp_geneid_param = "$work_dir/$species.geneid.param.tmp";
                        
                        #~ $param->writeParam($temp_geneid_param);
                        #~ my $my_command;

                        #~ #my $fh_geneid;
                        #~ #my $fname_geneid;
                        #~ print "\nbefore running my_command(s) geneid \n";
                        #~ my $fh_geneid    = File::Temp->new();
                        #~ my $fname_geneid = $fh_geneid->filename;
                        #~ print "\ntemp geneid file: $fname_geneid \n";
                        #~ $my_command =
#~ "./bin/geneid -GP $temp_geneid_param $gpfa > $fname_geneid";
                        #~ run($my_command);
                        #~ my $fh_gawk_out    = File::Temp->new();
                        #~ my $fname_gawk_out = $fh_gawk_out->filename;
                        #~ $my_command =
#~ "cat $fname_geneid | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' > $fname_gawk_out";
                        #~ run($my_command);
                        #~ $my_command =
#~ "egrep -wv 'exon'  $fname_gawk_out > $tmp_dir/Predictions.${newparam}.gff";
                        #~ run($my_command);

#~ #` ./bin/geneid -GP ${newparam}.temp $gpfa | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' > $tmp_dir/Predictions.${newparam}.gff`;
                        #~ print "\nafter running my_command(s) geneid \n";

                        #~ my @evaluation_output = split " ",
#~ ` ./bin/evaluation -sta $tmp_dir/Predictions.${newparam}.gff $gpgff | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $IoWF, $IeWF, $iAccCtx, $iMin, \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}'`;

                        #~ push( @evaluation_total, \@evaluation_output );

                    #~ }    #$iAccCtx for loop
                    #~ unless ( $iMin == $fMin ) {
                        #~ print STDERR "\nMinDist:  ";
                    #~ }
                    #~ $iAccCtx = $iAccCtxini;

                #~ }    #$iMin for loop
                #~ unless ( $IoWF == $FoWF ) {
                    #~ print STDERR "\noWF:  ";
                #~ }
                #~ $iAccCtx = $iAccCtxini;
                #~ $iMin    = $iMinini;

            #~ }    #if $branchswitch
            #~ elsif ( !$branchswitch ) 


sub sortevalbranch {

         $b->[9] <=> $a->[9]
      || $b->[12] <=> $a->[12]
      || $b->[6] <=> $a->[6]
      || $a->[13] <=> $b->[13]
      || $a->[14] <=> $b->[14]

}

   #~ elsif ($branchswitch) {

        #~ @sortedeval = sort sortevalbranch @$evalarray;

        #~ $bestIoWF = $sortedeval[0][0];    #0.2
        #~ $bestIeWF = $sortedeval[0][1];    #-3.5
        #~ $bestAcc  = $sortedeval[0][2];
        #~ $bestMin  = $sortedeval[0][3];

        #~ print STDERR "\nBest performance obtained using IoWF: "
          #~ . $sortedeval[0][0]
          #~ . " , IeWF: "
          #~ . $sortedeval[0][1]
          #~ . ", MinimalBranchDist: "
          #~ . $bestMin
          #~ . ", Optimal Branch Context: "
          #~ . $bestAcc . "\n";
        #~ print $fh_SOUT "best performance obtained using IoWF: "
          #~ . $sortedeval[0][0]
          #~ . " , IeWF: "
          #~ . $sortedeval[0][1]
          #~ . ", MinimalBranchDist: "
          #~ . $bestMin
          #~ . ", Optimal Branch Context: "
          #~ . $bestAcc . "\n";

        #~ #INITIALIZE ARRAY WITH EVALUATION PARAMETERS
        #~ @evaluationinit =
          #~ (
            #~ qw(oWF eWF AccCtx MinD SN SP CC SNe SPe SNSP SNg SPg SNSPg raME raWE)
          #~ );

        #~ print STDERR
#~ "(Sorted performance results (Three best performance estimates) for different values of oWF, eWF, AccCtx and MinD)\n\n"
          #~ . join( "\t", @evaluationinit ), "\n";
        #~ print $fh_SOUT
#~ "(Sorted performance results (best to worst) for different values of oWF, eWF, AccCtx and MinD)\n\n"
          #~ . join( "\t", @evaluationinit ), "\n";

        #~ foreach my $eval_ref (@sortedeval) {

            #~ print $fh_SOUT join( "\t", @$eval_ref ), "\n";

        #~ }

#~ ###THREE BEST PERFORMANCE RESULTS AFTER OPTIMIZATION
        #~ for ( my $i = 0 ; $i <= 2 ; $i++ ) {
            #~ print STDERR join( "\t", @{ $sortedeval[$i] } ), "\n";
        #~ }
#~ ############

        #~ ## BUILD BEST PERFORMANCE (OPTIMIZED) PARAMETER FILE (USE VALUES OBTAINED ABOVE)

        #~ my $param = Geneid::Param->new();
        #~ $param->readParam("$newparam");

        #~ for ( my $i = 0 ; $i < $param->numIsocores ; $i++ ) {
            #~ if (
                #~ !defined @{ $param->isocores }[$i]->Exon_weights(
                    #~ [ $bestIeWF, $bestIeWF, $bestIeWF, $bestIeWF ]
                #~ )
              #~ )
            #~ {
                #~ croak "error in setting exon weights\n";
            #~ }
            #~ if ( !defined @{ $param->isocores }[$i]
                #~ ->Exon_factor( [ $bestIoWF, $bestIoWF, $bestIoWF, $bestIoWF ] )
              #~ )
            #~ {
#~ #   if (!defined @{$param->isocores}[$i]->Exon_factor([0.4,$bestIoWF,$bestIoWF,0.4])) {
                #~ croak "error in setting exon weights\n";
            #~ }
            #~ if (
                #~ !defined @{ $param->isocores }[$i]->Site_factor(
                    #~ [
                        #~ 1 - $bestIoWF,
                        #~ 1 - $bestIoWF,
                        #~ 1 - $bestIoWF,
                        #~ 1 - $bestIoWF
                    #~ ]
                #~ )
              #~ )
            #~ {
#~ #  if (!defined @{$param->isocores}[$i]->Site_factor([0.55,1-$bestIoWF,1-$bestIoWF,0.55])) {
                #~ croak "error in setting exon weights\n";
            #~ }
            #~ if (
                #~ !defined @{ $param->isocores }[$i]->set_profile(
                    #~ 'Branch_point_profile', $prof_len_bra, $fxdbraoffset,
                    #~ -50, 0, 0, 1, $bestAcc, $bestMin, 0, 0, $branchmatrix
                #~ )
              #~ )
            #~ {
                #~ croak "error in setting profile\n";
            #~ }
        #~ }

        #~ #write new parameter file (optimized)
        #~ $param->writeParam("$species.geneid.optimized.param");

        #~ print STDERR
#~ "\nNew optimized parameter file named: $species.geneid.optimized.param \n";
        #~ print $fh_SOUT
#~ "\nNew optimized parameter file named: $species.geneid.optimized.param \n";

        #~ return [ $bestIeWF, $bestIoWF, $bestAcc, $bestMin, \@evaluationinit ];

    #~ }    #if branch switch
    
