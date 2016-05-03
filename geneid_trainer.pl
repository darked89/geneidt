#!/usr/bin/perl

## if run under perlbrew, use i.e.:
#!/$HOME/perl5/perlbrew/perls/perl-5.10.1/bin/perl

## checks & debugs modules
use Modern::Perl;
use English '-no_match_vars';

use strict;
use warnings;
use autodie;
use diagnostics -verbose;
use sigtrap qw(stack-trace old-interface-signals);

use Carp qw(carp cluck croak confess);
use Carp::Always;
use Carp::Assert qw(assert);

#use Data::Dumper::Perltidy;
use Data::Dumper;

## common modules, used in part at the moment
#use Env qw(PATH, PERL5LIB);
use Cwd;
use Getopt::Long;

use File::Path;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;

use IPC::System::Simple qw(run system capture EXIT_ANY);
use Readonly;

use feature 'say';
use Benchmark qw(:all);
use Devel::Size qw(size total_size);

## geneid_trained modules
use Geneid::Param;
use Geneid::Isocore;
use Geneid::geneid;
use Geneid::geneidCEGMA;

## experimental
#~ use Inline::Python;

## MAIN VARIABLES
my $PROGRAM      = "geneid_trainer";
my $VERSION      = "2016.04.26";
my $PROGRAM_HOME = getcwd;

my $exec_path = "$PROGRAM_HOME/bin/";

local $ENV;
$ENV{'PATH'} = $exec_path . ":" . $ENV{'PATH'};

my $genetic_code = "./etc/genetic.code";

#print STDERR "geneticcode: $genetic_code\n";

## no need to run anything if this fails
check_external_progs();

## Move parts necessary for getting comand line args here
my $species      = "";
my $input_gff_fn = "";
my $input_fas_fn = "";
my $sout         = "-";
my $fix_fasta    = 0;

my $branchp = 0;

#my $reduced     = 0;
my $interactive = 0;

#~ my $tenfold     = 0;
#~ my $input_gff_fn2ps      = 0;
#~ my $pout        = "-";

my $standard_run = 1;    # get rid of !$reduced

## Get arguments (command line)
GetOptions(
    'species:s'       => \$species,
    'gff:s'           => \$input_gff_fn,
    'fasta:s'         => \$input_fas_fn,
    'sout|statsout:s' => \$sout,
    'fix_fasta'       => \$fix_fasta

      #'branch'          => \$branchp,
      #'reduced|red'     => \$reduced,

      #'path|binpath:s'  => \$path,
      #'interactive|i' => \$interactive,
      #'tenfold'       => \$tenfold,
      #'gff2ps'        => \$input_gff_fn2ps
);
my $usage =
"Usage: $PROGRAM_NAME -species H.sapiens -gff gff_name -fasta fasta_name -sout statsfile_out -fix_fasta 1";

#~ -branch -reduced -path <executables_path>\n";

if ( !( $species && $input_gff_fn && $input_fas_fn && $sout && $fix_fasta ) ) {
    print STDERR $usage and exit;
}

#  unless ( $species && $input_gff_fn && $input_fas_fn && $sout && $fix_fasta );
## EXAMPLE COMMAND LINE: ./geneidTRAINer1_2TA.pl -species S.cerevisiae -gff S_cerevisiae4training.gff -fastas yeast_genome.fa -sout stats.txt -branch -reduced
## Get arguments (command line) END

## end of getting ARGS
## set CONSTANTS
## Constant values. modify if needed
Readonly::Scalar my $pwm_cutoff     => -7;
Readonly::Scalar my $bases_offset   => 30;    #bases in front/after? a feature
Readonly::Scalar my $train_fraction => 0.8;   #fraction of seq used for training
Readonly::Scalar my $train_loci_cutoff         => 500;
Readonly::Scalar my $train_sites_cutoff        => 1400;
Readonly::Scalar my $train_sites_cutoff_alt    => 1200;  # changed in some part?
Readonly::Scalar my $train_sites_markov_cutoff => 5500;
Readonly::Scalar my $backgrnd_kmer_size        => 60;
Readonly::Scalar my $backgrnd_kmer_num => 100_000;

#~ ( $totalcodingbases > 400000 && $totalnoncodingbases > 100000 )
Readonly::Scalar my $coding_bp_limit_A => 400_000;
Readonly::Scalar my $coding_bp_limit_B => 375_000;

#~ || ( $totalcodingbases > 375000 && $totalnoncodingbases > 150000 )
Readonly::Scalar my $non_coding_bp_limit_A => 100_000;
Readonly::Scalar my $non_coding_bp_limit_B => 150_000;

## need to explain or just incorporate Readonly::Scalar
## BUG ??? is it ever used?
##  ## Q_FRANCISCO: condition below?
#~ || (   $totalnoncodingbases > 35000
#~ && $totalcodingbases > ( 25 * $totalnoncodingbases ) )

my %profile_params = (
## EXON WEIGHT PARAMETER
    IeWF => -4.5,
    deWF => 0.5,
    FeWF => -2.5,
## EXON/OLIGO FACTOR PARAMETER
    IoWF => 0.25,
    doWF => 0.05,
    FoWF => 0.50,
## Minimum Branch Profile Distance
    iMin => 7,
    dMin => 2,
    fMin => 9,
## ACCEPTOR CONTEXT
    iAccCtx => 40,
    dAccCtx => 10,
    fAccCtx => 70,
);

## end set CONSTANTS

## hanging to /tmp for faster exec on clusters
## TODO: random strin in dir name to avoid conflicts with other ppl runnig the script
my $work_dir = "/tmp/workdir_00_gtrain/";

my $tmp_dir     = "$work_dir/temp_00/";
my $stats_dir   = "$work_dir/stats/";
my $sites_dir   = "$work_dir/sites/";
my $plots_dir   = "$work_dir/plots/";
my $introns_dir = "$work_dir/introns/";
my $fastas_dir  = "$work_dir/fastas/";
my $backgrd_dir = "$work_dir/backgrd/";
my $cds_dir     = "$work_dir/cds/";
my $geneid_dir  = "$work_dir/geneid/";
my $results_dir = "$work_dir/results/";

my @data_dirs = (
    $work_dir,  $tmp_dir,    $fastas_dir,  $stats_dir,
    $sites_dir, $plots_dir,  $introns_dir, $backgrd_dir,
    $cds_dir,   $geneid_dir, $results_dir
);

create_data_dirs(@data_dirs);

## PROGRAM SPECIFIC VARIABLES (unordered...)

my $fh_SOUT;
my $last_bench_time;    #for benchmarking parts of the script
my $no_dots_gff_fn = "$work_dir/input_gff_no_dots.gff";
say "\nXXX $no_dots_gff_fn \n";

my $temp_GENOMEX_tblcaps = "";
my $backgrnd_kmers_fn    = "";
## my $tblseq      = "";

my @evaluation = ();
## flow control / option variables
my $run_jacknife_flag   = 0;
my $use_allseqs_flag    = 0;
my $use_branch_flag     = 0;
my $run_optimizeX_flag  = 0;    ## UNUSED
my $run_contig_opt_flag = 0;
## my $ext_flg          = 0;

my $tmp_4train_gff      = "";
my $tmp_4eval_gff       = "";
my $tmp_locus_id_X      = "";
my $tmp_locus_id_X_new  = "";
my $tmp_locus_id_eval_X = "";
my $seqs_used_XX        = "";

## Golden Path stuff
my $gp_evalcontig_fa      = "";
my $gp_evalcontig_gff     = "";
my $gp_evalcontig_len_int = "";
my $gp_evalcontig_tbl     = "";
my $gp_eval_fa_X          = "";
my $gp_eval_gff_X         = "";
my $gp_eval_len_intX      = "";
my $gp_eval_tbl_X         = "";

my $gp_traincontig_fa      = "";
my $gp_traincontig_gff     = "";
my $gp_traincontig_len_int = "";
my $gp_traincontig_tbl     = "";
my $gp_train_fa_X          = "";
my $gp_train_gff_X         = "";
my $gp_train_len_int_X     = "";
my $gp_train_tbl_X         = "";

## Weights for profiles??
my $best_IeWF = "";
my $best_IoWF = "";
my $best_Acc  = "";
my $best_Min  = "";

my $subprofile_acceptors = "";
my $subprofile_donors    = "";
my $subprofile_ATGs      = "";

my $X_geneid_sorted_gff = "";
my $seqs_eval_gff_X     = "";
my $inframe_X           = 0;
my $inframe_X_eval      = 0;
my $locus_id            = "";
my $locus_id_new        = "";
my $intron_long_int     = 0;
my $intron_short_int    = 0;
my $intergenic_max      = "";
my $intergenic_min      = "";

my $out_acceptor_tbl       = "";
my $cds_all_nozero_tbl     = "";
my $cds_eval_nonzero_tbl   = "";
my $out_donor_tbl          = "";
my $out_gff_X              = "";
my $out_eval_gff_X         = "";
my $out_intron_X           = "";
my $out_intron_eval_X      = "";
my $out_locus_id_X         = "";
my $out_locus_id_X_eval    = "";
my $out_ATGxs_tbl          = "";
my $seqs_4evaluation_listX = "";
my $seqs_4training_listX   = "";

my $ATGx_tbl                   = "";
my $tmp_X_geneid_sorted_gff    = "";
my $train_2cols_seq_locusid_fn = "";
my $eval_2cols_seq_locusid_fn  = "";
my $templist_train             = "";

my $tot_noncanon_accept_intX = "";
my $tot_noncanon_donors_intX = "";
my $tot_noncanon_ATGx        = "";

my $total_seqs             = "";
my $tot_seqs4training_intX = "";

#############################################################
## INITIAL CHECKS
#############################################################
## TODO  fasta / gff file accessible?
sub validate_input_fasta {
    return 1;
}

sub validate_input_gff {
    return 1;
}

sub select_test_eval {
    return 1;
}

sub extract_CDS {
    return 1;
}

sub extract_introns {
    return 1;
}

sub extract_donors {
    return 1;
}

sub extract_acceptors {
    return 1;
}

sub extract_ATGx {
    return 1;
}

## TODO 2. limits:
## 2a. >= 500 genes in gff
## PYTHON INLINE START
# use Inline Python => << 'PYEND';
# from pygeneid import check_fasta

# PYEND

# my $headers_fasta_seq;
# $headers_fasta_seq  = check_fasta($input_fas_fn);

# print STDERR "\n PYTHON The user has provided $headers_fasta_seq  genomic sequences\n";

##die;
## PYTHON INLINE END
#############################################################
## Common tasks
#############################################################
## sanity check
#######################################################
## CREATE FASTAS CDS; INTRON, SITES DIRs WITHIN PATH (ONLY FIRST TIME)

## CREATE BLANK PARAMETER FILE############
my $param    = Geneid::Param->new($species);
my $newparam = "$work_dir/$species.geneid.param";

#set isochores to 1
$param->numIsocores(1);
$param->isocores( [ Geneid::Isocore->new() ] );

## END CREATING A PARAMETER FILE REGARDLESS OF WHETHER THE TRAINING IS COMPLETE OR REDUCED

normal_run();

sub normal_run {
    my $fh_FOUT;
    my $my_command;
    my $old_option = 0;

## Convert fasta to tabular format
## Fasta process to sub later
    my $t0 = Benchmark->new;

    if ($old_option) {
        print "not yet\n";
    }
    print STDERR
      "\nConverting genomics fasta file ($input_fas_fn) to tabular format\n";

    my $genomic_temp_tbl = $work_dir . $species . ".genomic.tmp.tbl";
    fasta_2_tbl( $input_fas_fn, $genomic_temp_tbl );
    run("sort -o $genomic_temp_tbl $genomic_temp_tbl");

    print STDERR "actg to ACTG conversion of input fasta \n";
    my $tblcaps = "";

    open(
        my $fh_LOCID,
        "-|",
"gawk '{gsub(/_/,\"\",\$1);gsub(/\\./,\"\",\$1);print \$1, toupper(\$2)}' $genomic_temp_tbl "
    ) or croak "Failed here";
    while (<$fh_LOCID>) {
        $tblcaps .= $_;
    }
    close $fh_LOCID;

    chomp $tblcaps;
    $temp_GENOMEX_tblcaps = $work_dir . $species . ".genomic.tbl";
    open( my $fh_FOUT_caps, ">", "$temp_GENOMEX_tblcaps" )
      or croak "Failed here";
    print {$fh_FOUT_caps} "$tblcaps";
    close $fh_FOUT_caps;

## place genomic sequences in "fastas_$species" directory
    print STDERR "move genomic sequences into \"$fastas_dir\" directory\n";
    print STDERR "(also transfer genomic fasta length info)\n\n";
## do not create fastas in diretory if they are already created and their number corresponds to the number of sequences in thr array
## CONVERT GENOMICS FASTA TO MULTI FASTA AND PLACE THEM IN APPROPRIATE DIRECTORY

    print STDERR
"Convert $temp_GENOMEX_tblcaps to multiple genomic fastas and place them in $fastas_dir:\n";

    TblToFastaFile( $fastas_dir, $temp_GENOMEX_tblcaps );
    print STDERR
"\n\nConversion of $temp_GENOMEX_tblcaps to multiple genomic fastas completed..\n\nAdd fasta sequence length information to same directory\n\n";
    write_sizes_from_tbl_fn($temp_GENOMEX_tblcaps);

#################################################

## get locus_id file only first time pipeline is run for a given species #ALL GENE MODELS
    ##  ## Q_FRANCISCO: can we assume some sane GFF/GFT format as an input?  Why _ and "."?

    if ($old_option) {
        print STDERR
          "\nEliminate undesirable (_ and .) characters from $input_gff_fn\n";

        my $filtergff = "";

        $my_command =
"gawk '{OFS=\"\\t\"}{gsub(/\\./,\"\",\$1);gsub(/\\./,\"\",\$9);gsub(/_/,\"\",\$0);print}' $input_gff_fn";
        $filtergff = capture($my_command);

        open( my $fh_FOUT, ">", "$no_dots_gff_fn" ) or croak "Failed here";
        print $fh_FOUT "$filtergff";
        close $fh_FOUT;

    }
    else {
        $no_dots_gff_fn = $input_gff_fn;
    }
    print STDERR "\nObtain locus_id (list of genomic sequences / genes)\n";

    $my_command = "gawk '{print \$1,\$9}' $no_dots_gff_fn | sort | uniq ";
    $locus_id   = capture($my_command);

    # say "\n TTT got here TTT\n";
    $tmp_locus_id_X = $work_dir . $species . "_locus_id";
    open( $fh_FOUT, ">", "$tmp_locus_id_X" ) or croak "Failed here";
    print $fh_FOUT "$locus_id";
    close $fh_FOUT;

## number of gene models TOTAL
    $my_command = " gawk '{print \$2}' $tmp_locus_id_X | sort | uniq | wc -l";
    $total_seqs = capture($my_command);

#~ $total_seqs =
#~ ` gawk '{print \$2}' $tmp_locus_id_X | sort | uniq | wc | gawk '{print \$2}' `;
    chomp $total_seqs;
## number of genomic sequences TOTAL
    $my_command = "gawk '{print \$1}' $tmp_locus_id_X | sort | uniq | wc -l";
    my $total_genomic = capture($my_command);

#~ $total_genomic =
#~ ` gawk '{print \$1}' $tmp_locus_id_X | sort | uniq | wc | gawk '{print \$1}' `;
    chomp $total_genomic;

    print STDERR
"\nThe gff file ($no_dots_gff_fn) contains a total of $total_genomic genomic sequences and $total_seqs gene models\n";

## get a list of genes TOTAL
    print STDERR "\nObtain list of all genes\n\n";
    my $list_seqs = "";

    $my_command = "gawk '{print \$9}' $no_dots_gff_fn | sort | uniq ";
    $list_seqs  = capture($my_command);

    my $templist = $work_dir . $species . "_list_all_seqs";
    open( $fh_FOUT, ">", "$templist" ) or croak "Failed here";
    print $fh_FOUT "$list_seqs";
    close $fh_FOUT;

    if ( $total_seqs >= $train_loci_cutoff ) {

        $tot_seqs4training_intX = int( $train_fraction * $total_seqs );

        print STDERR
"\nA subset of $tot_seqs4training_intX sequences (randomly chosen from the $total_seqs gene models) was used for training\n";
        ## DEBUG KEEP !!! shuf => random select
        ## head -$tot_seqs4training_intX just the first ones
#my $my_command =           "shuf --head-count=$tot_seqs4training_intX $tmp_locus_id_X | sort | uniq";

        my $my_command =
          "head --lines=$tot_seqs4training_intX $tmp_locus_id_X | sort | uniq";

        $locus_id_new = capture($my_command);

        $tmp_locus_id_X_new =
          $work_dir . $species . "_locus_id_training_setaside80";
        open( my $fh_FOUT, ">", "$tmp_locus_id_X_new" ) or croak "Failed here";
        print $fh_FOUT "$locus_id_new";
        close $fh_FOUT;

## ASSUMING USER SELECTED TO SET ASIDE SEQUENCES FOR EVALUATION (20%)
        $my_command =
          "gawk '{print \$2}' $tmp_locus_id_X_new | sort | uniq | wc -l";
        $seqs_used_XX = capture($my_command);
        chomp $seqs_used_XX;
        my $t1 = Benchmark->new;
        my $td = timediff( $t1, $t0 );
        print "\nTTT the code took t0->t1:", timestr($td), "\n";
        $last_bench_time = $t1;

###################
## gff for training subset
####################
        my $gff4training = "";

        print STDERR
"\nThe new training gff file includes $seqs_used_XX gene models (80% of total seqs)\n";
        ## ??? BUG ???
        $my_command =
"gawk '{print \$2\"\$\"}' $tmp_locus_id_X_new | sort | uniq | egrep -wf - $no_dots_gff_fn";
        $gff4training = capture($my_command);

        $tmp_4train_gff = $work_dir . $species . ".gff_training_setaside80";
        open( $fh_FOUT, ">", "$tmp_4train_gff" ) or croak "Failed here";
        print $fh_FOUT "$gff4training";
        close $fh_FOUT;

        print STDERR "\nObtain list of training genes\n\n";

        my $list_seqs_train = "";

        $my_command      = "gawk '{print \$9}' $tmp_4train_gff | sort | uniq ";
        $list_seqs_train = capture($my_command);

        $templist_train = $work_dir . $species . "_list_train_seqs_setaside80";
        open( $fh_FOUT, ">", "$templist_train" ) or croak "Failed here";
        print $fh_FOUT "$list_seqs_train";
        close $fh_FOUT;

#########################
## new locus_id for evaluation test set
#########################
        my $locusideval = "";

        $my_command =
"gawk '{print \$0\"\$\"}' $templist_train | egrep -vwf - $tmp_locus_id_X";
        $locusideval = capture($my_command);
        chomp $locusideval;

        $tmp_locus_id_eval_X =
          $work_dir . $species . "_locus_id_evaluation_setaside20";
        open( $fh_FOUT, ">", "$tmp_locus_id_eval_X" ) or croak "Failed here";
        print $fh_FOUT "$locusideval";
        close $fh_FOUT;

######################
## gff for evaluation test set
#########################

        $my_command =
"gawk '{print \$2\"\$\"}' $tmp_locus_id_eval_X | sort | uniq | egrep -wf - $no_dots_gff_fn | gawk '{ print \$9}' | sort | uniq | wc -l";
        ## ??? BUG this is a number...
        $seqs_eval_gff_X = capture($my_command);

        #chomp $input_gff_fnseqseval;

        print STDERR
"The evaluation gff file includes $seqs_eval_gff_X gene models (20% of total seqs)\n\n";

        $my_command =
"gawk '{print \$2\"\$\"}' $tmp_locus_id_eval_X | sort | uniq | egrep -wf - $no_dots_gff_fn ";
        my $gff4evaluation = capture($my_command);

        $tmp_4eval_gff = $work_dir . $species . ".gff_evaluation_setaside20";
        open( $fh_FOUT, ">", "$tmp_4eval_gff" ) or croak "Failed here";
        print $fh_FOUT "$gff4evaluation";
        close $fh_FOUT;

    }    # seqs > 500
####LOOP IF WE HAVE FEWER THAN 500 SEQUENCES

    else {    # seqs < $train_loci_cutoff
        ## BUG we do not do jacknife anyway here
        croak "we do not have >= $train_loci_cutoff sequences, quitting now";
    }    # seqs < 500

    if ( !$use_allseqs_flag ) {    ##SET SEQS FOR EVAL AND TRAINING (SUBSETS)
        print STDERR
          "\nConvert general gff2 to geneid-gff format  NOT_USE_ALL_SEQS \n\n";
        ### XXX function name
    }
## Convert general gff2 to geneid gff format
## extract and check cds and intron sequences. Remove inframe stops and check all seqs start with ATG and end with STOP
## TRAIN
    #~ my $t1 = Benchmark->new;
    #~ my $td = timediff($t1, $t0);
    #~ print "\nTTT the code took t0->t1:",timestr($td),"\n";
    #~ $last_bench_time = $t1;

## XXXXXX extracted lines start
    {

## BUG => there is no need to convert gff to geneid format each time we run

        $train_2cols_seq_locusid_fn =
          convert_GFF_2_geneidGFF( $tmp_4train_gff, $species, ".train" );
        print STDERR
          "L575 : $train_2cols_seq_locusid_fn \t $tmp_locus_id_X_new \n";

        (
            $cds_all_nozero_tbl, $out_intron_X, $out_locus_id_X, $out_gff_X,
            $inframe_X
          )
          = @{
            extractCDSINTRON( $train_2cols_seq_locusid_fn, $tmp_locus_id_X_new,
                ".train" )
          };

#  print STDERR " OUTSIDE EXTRACTCDSINTRON outgff: $out_gff_X\noutlocus_id: $out_locus_id_X\n";
## TRAIN
## EVAL
        $eval_2cols_seq_locusid_fn =
          convert_GFF_2_geneidGFF( $tmp_4eval_gff, $species, ".eval" );

        print STDERR
"L588 tmp_locus_id_X_new:  $eval_2cols_seq_locusid_fn \t $tmp_locus_id_eval_X,\n";
        (
            $cds_eval_nonzero_tbl, $out_intron_eval_X, $out_locus_id_X_eval,
            $out_eval_gff_X,       $inframe_X_eval
          )
          = @{
            extractCDSINTRON( $eval_2cols_seq_locusid_fn, $tmp_locus_id_eval_X,
                ".eval" )
          };
###EVAL

    }

###########

## extract and check splice sites and start codon. Use only canonical info #IN SEQUENCES USED IN TRAINING
    print STDERR "L623 :  $out_gff_X \t $out_locus_id_X  \n";
    (
        $out_donor_tbl,    $tot_noncanon_donors_intX,
        $out_acceptor_tbl, $tot_noncanon_accept_intX,
        $out_ATGxs_tbl,    $tot_noncanon_ATGx
    ) = @{ extractprocessSITES( $out_gff_X, $out_locus_id_X ) };

## prepare sequences for optimization of newly developed parameter file (TRAIN)

    print STDERR
"\nConvert gff to gp (golden-path-like)format (artificial contig - concatenated sequences - approx. 800 nt between sequences)\n";

    (
        $gp_traincontig_gff, $gp_traincontig_fa,
        $gp_traincontig_tbl, $gp_traincontig_len_int
    ) = @{ processSequences4Optimization( $out_gff_X, ".train", 1 ) };

    print STDERR
"\nConvert gff to gp (golden-path-like)format (training set for later optimization -400-nt flanked sequences)\n";
    ( $gp_train_gff_X, $gp_train_fa_X, $gp_train_tbl_X, $gp_train_len_int_X ) =
      @{ processSequences4Optimization( $out_gff_X, ".train", 0 ) };
    print STDERR "$gp_train_gff_X";

########################################

## NOT USING ALL SEQS FOR TRAINING/EVALUATION ALSO PROCESS EVAL SEQS
    if ( !$use_allseqs_flag ) {
        print STDERR
"\n NNN NOT USING ALL SEQS FOR TRAINING/EVALUATION ALSO PROCESS EVAL SEQS\nn";

## prepare test set for evaluation of newly developed parameter file (EVAL)

        print STDERR
"\nConvert gff to gp (golden-path-like)format (400-nt flanking)(test set for evaluation of new parameter file)\n";
        ( $gp_eval_gff_X, $gp_eval_fa_X, $gp_eval_tbl_X, $gp_eval_len_intX ) =
          @{ processSequences4Optimization( $out_eval_gff_X, ".eval", 0 ) };
        print STDERR "DONE\n";

        print STDERR
"\nConvert gff to gp (golden-path-like)format (test set for evaluation of new parameter file - (artificial contig - concatenated sequences - approx. 800 nt between sequences)\n";
        (
            $gp_evalcontig_gff, $gp_evalcontig_fa,
            $gp_evalcontig_tbl, $gp_evalcontig_len_int
        ) = @{ processSequences4Optimization( $out_eval_gff_X, ".eval", 1 ) };
        print STDERR "DONE\n";

    }

    my $t3 = Benchmark->new;
    my $td = timediff( $t3, $last_bench_time );
    print "\nTTT the code took t2->t3:", timestr($td), "\n";
    $last_bench_time = $t3;

## XXXXXX extracted lines start

### EVERYTHING BELOW ALWAYS EXECUTED EVEN ON SHORT VERSION OF THE PIPELINE (REDUCED)

## GET BACKGROUND SEQUENCES

    print
"Obtaining $backgrnd_kmer_num background sequences of $backgrnd_kmer_size nucleotides each for estimating background frequencies of nucleotides\n";

    $backgrnd_kmers_fn = $work_dir . $species . "_background.tbl";
    get_background_kmers(
        $backgrnd_kmer_size, $input_fas_fn, $temp_GENOMEX_tblcaps,
        $backgrnd_kmer_num,  $backgrnd_kmers_fn
    );

    my (
        $donor_start,  $donor_end,  $acceptor_start,
        $acceptor_end, $ATGx_start, $ATGx_end
    ) = compute_sites_pictogram();

    #~ ## DEBUG MEM

    #~ {
    #~ #no scrict 'refs';
    #~ my $size = 0;
    #~ for my $var (keys %{'main::'}) {
    #~ $size = size("A string");
    #~ print "$var = $size \n";
    #~ }
    #~ }
    #~ ## DEBUG MEM END

## DERIVE INITIAL/TRANSITION MARKOV MODEL

    my (
        $markov_mod_ini,  $markov_mod_trans, $total_coding,
        $total_noncoding, $markov_model
    ) = @{ deriveCodingPotential( $cds_all_nozero_tbl, $out_intron_X ) };

    #add markov matrices to the parameter file
    if ( !defined @{ $param->isocores }[0]->Markov_order($markov_model) ) {
        croak "error in setting Markov_order\n";
    }
    if ( !defined @{ $param->isocores }[0]
        ->Markov_Initial_probability_matrix($markov_mod_ini) )
    {
        croak "error in setting Markov_Initial_probability_matrix\n";
    }
    if ( !defined @{ $param->isocores }[0]
        ->Markov_Transition_probability_matrix($markov_mod_trans) )
    {
        croak "error in setting Markov_Transition_probability_matrix\n";
    }
######################################

    ## BUG do not remove, misleading name of a function: WriteStatsFile
    ( $intron_short_int, $intron_long_int, $intergenic_min, $intergenic_max ) =
      calculate_stats(
        $species,                  $sout,
        $out_intron_X,             $cds_all_nozero_tbl,
        $out_gff_X,                $inframe_X,
        $inframe_X_eval,           $seqs_used_XX,
        $tot_noncanon_donors_intX, $tot_noncanon_accept_intX,
        $tot_noncanon_ATGx,        $markov_model,
        $total_coding,             $total_noncoding,
        $donor_start,              $donor_end,
        $acceptor_start,           $acceptor_end,
        $ATGx_start,               $ATGx_end,
        0,                         0,
        0,                         $use_allseqs_flag
      );

    print STDERR
"\nshortest intron: $intron_short_int\nlongest intron: $intron_long_int\nminimum intergenic: $intergenic_min\nmaximum intergenic: $intergenic_max\n";

##################################################
## WRITE PRELIMINARY NON-OPTIMIZED PARAMETER FILE

    my $newparam = "$work_dir/$species.geneid.param";
    $param->writeParam($newparam);

### if reduced training (non-default) do not calculate any of the above ALL OF THE ABOVE MUST BE RUN ONLY FIRST TIME GENEID IS TRAINED FOR A GIVEN SPECIES
###EVERYTHING BELOW WILL BE RUN EVERYTIME THE TRAINING PIPELINE IS RUN WHETHER "REDUCED" OR "FULL"

################################
## OPTIMIZE PARAMETER FILE
################################

    print STDERR "\nOptimizing new parameter file\n\n";

## BUG we run non interactive here
## simplification
    my $opttype = "";
    $opttype             = "contig";
    $run_contig_opt_flag = 1;
    $run_jacknife_flag   = 0;

## BUG settings need to be set forward. Also these are numbers.

## OPTIMIZATION FUNCTION NO BRANCH
    my $array_ref = "";

## EXON WEIGHT PARAMETER
    my $IeWF = -4.5;
    my $deWF = 0.5;
    my $FeWF = -2.5;
## EXON/OLIGO FACTOR PARAMETER
    my $IoWF = 0.25;
    my $doWF = 0.05;
    my $FoWF = 0.50;
##Minimum Branch Profile Distance
    my $iMin = 7;
    my $dMin = 2;
    my $fMin = 9;
##ACCEPTOR CONTEXT
    my $iAccCtx = 40;
    my $dAccCtx = 10;
    my $fAccCtx = 70;

## OPTIMIZATION FUNCTIONS

    if ( !$run_contig_opt_flag ) {
        print STDERR "\n DEBUG: NOT CONTIG OPT\n";
        croak "bad bug!\n";
        ;
    }    #end if
         #~ if ( !$run_contig_opt_flag ) {
         #~ print STDERR "\n DEBUG: NOT CONTIG OPT\n";
         #~ @evaluation = @{
         #~ OptimizeParameter( $gp_train_fa_X, $gp_train_gff_X, $newparam, 0,
         #~ 0, 0, 0,
         #~ $IeWF, $deWF, $FeWF, $IoWF, $doWF, $FoWF, 0, 0, 0, 0, 0, 0 )
         #~ };

    #~ ( $best_IeWF, $best_IoWF, $best_Acc, $best_Min, $array_ref ) = @{
    #~ BuildOptimizedParameterFile( \@evaluation, $use_branch_flag, 0, 0,
    #~ 0 )
    #~ };

    #~ }
    elsif ($run_contig_opt_flag) {
        print STDERR "\n DEBUG: CONTIG OPT\n";
        @evaluation = @{
            OptimizeParameter( $gp_traincontig_fa, $gp_traincontig_gff,
                $newparam, 0, 0, 0, 0, $IeWF, $deWF, $FeWF, $IoWF, $doWF,
                $FoWF, 0, 0, 0, 0, 0, 0 )
        };

        ( $best_IeWF, $best_IoWF, $best_Acc, $best_Min, $array_ref ) = @{
            BuildOptimizedParameterFile( \@evaluation, $use_branch_flag, 0, 0,
                0 )
        };

    }    #elsif end

    my @evaluationinit = @{$array_ref};
    my @evaluationtest = ();

############
## EVALUATE PERFORMANCE OF NEW PARAMETER FILE ON TEST SET (IF PRESENT)
############

    my $paramopt = "$species.geneid.optimized.param";

    if ( !$use_allseqs_flag ) {
        my $fh_SOUT;
        open( $fh_SOUT, ">", "$work_dir/$species.use_NOT_allseqs.log" );

     #print STDERR "CHECK EVALUATE: $gp_eval_fa_X, $gp_eval_gff_X, $paramopt\n";

        if ( !$run_contig_opt_flag ) {

            @evaluationtest = @{
                EvaluateParameter( $gp_eval_fa_X, $gp_eval_gff_X, $paramopt,
                    $IoWF, $IeWF )
            };

        }    #end if !$run_contig_opt_flag

        elsif ($run_contig_opt_flag) {

            @evaluationtest = @{
                EvaluateParameter(
                    $gp_evalcontig_fa, $gp_evalcontig_gff, $paramopt,
                    $IoWF,             $IeWF
                )
            };
        }

        if ( !$use_branch_flag ) {

            print STDERR
              "\nPerformance of new optimized parameter file on test set:\n\n"
              . join( "\t", @evaluationinit[ 2 .. $#evaluationinit ] ), "\n";

        }
        elsif ($use_branch_flag) {

            print STDERR
              "\nPerformance of new optimized parameter file on test set:\n\n"
              . join( "\t", @evaluationinit[ 4 .. $#evaluationinit ] ), "\n";

        }

        print STDERR join( "\t", @evaluationtest ), "\n\n";

        print $fh_SOUT join( "\t", @evaluationtest ), "\n\n";
        close $fh_SOUT;
    }    # if NOT using all seqs for training

    return 1;
} ## end normal run

#######################################
## END OF MAIN PORTION OF SCRIPT
#######################################

sub extractCDSINTRON {

    my ( $no_dots_gff_fn, $locus_id, $type ) = @_;

    # #####extract CDS and INTRON SEQUENCES
    #my
    print STDERR "\nEXTRACT CDS and INTRON SEQUENCES from $type set..\n\n";
    open( my $fh_LOCUS, "<", "$locus_id" ) or croak "Failed here";
    print STDERR "$locus_id and $no_dots_gff_fn\n";
    my $count = 0;
    while (<$fh_LOCUS>) {
        my ( $genomic_id, $gene_id ) = split;
        run(" egrep -w '$gene_id\$' $no_dots_gff_fn > $tmp_dir/$gene_id.gff");
        ## POTENTIAL BUG, split commands below

        #~ my $fh_ssgff_A    = File::Temp->new();
        #~ my $fname_ssgff_A = $fh_ssgff_A->filename;
        ## CDS
        my $fname_ssgff_A = "$cds_dir/${species}_ssgff_A.tmp.fa";
        my $my_command =
"./bin/ssgff -cE $fastas_dir/$genomic_id $tmp_dir/$gene_id.gff >  $fname_ssgff_A";
        run($my_command);
        $my_command =
"cat $fname_ssgff_A | sed -e 's/:/_/' -e 's/ CDS//' >> $cds_dir/${species}${type}.cds.fa ";
        run($my_command);

#` ./bin/ssgff -cE $work_dir/fastas_$species/$genomic_id $tmp_dir/$gene_id.gff | sed -e 's/:/_/' -e 's/ CDS//' >> $work_dir/cds/${species}${type}.cds.fa `;

        ## INTRONS
        #~ my $fh_ssgff_B    = File::Temp->new();
        #~ my $fname_ssgff_B = $fh_ssgff_B->filename;
        my $fname_ssgff_B = "$introns_dir/${species}_ssgff_B.tmp.fa";
        $my_command =
"./bin/ssgff -iE $fastas_dir/$genomic_id $tmp_dir/$gene_id.gff > $fname_ssgff_B";

        #say "\n$my_command\n";
        run($my_command);
        $my_command =
"cat $fname_ssgff_B | sed -e 's/:/_/' -e 's/ Intron.*//' >> $introns_dir/${species}${type}.intron.fa";

        #say "\n$my_command\n";
        run($my_command);

#` ./bin/ssgff -iE $work_dir/fastas_$species/$genomic_id $tmp_dir/$gene_id.gff | sed -e 's/:/_/' -e 's/ Intron.*//' >> $work_dir/intron/${species}${type}.intron.fa `;
        $count++;
        print STDERR "$count ..";
    }
    close $fh_LOCUS;

    print STDERR "DONE\n";

    # #####tabulate CDS and INTRON SEQUENCES

    print STDERR
"\nCreate tabular format of CDS and INTRON sequences for $type sequences\n";

## CDS
    my $cds_tmp_fa = $cds_dir . ${species} . "$type" . ".cds.fa";
    print STDERR "$cds_tmp_fa\n\n";
    my $cds_temp_tbl = $cds_dir . ${species} . "$type" . ".cds.tbl";
    fasta_2_tbl( $cds_tmp_fa, $cds_temp_tbl );
    print STDERR "cds tabular file created for $type sequences \n";

    # ##INTRON
    my $intron_tmp_fa  = $introns_dir . ${species} . "$type" . ".intron.fa";
    my $intron_tmp_tbl = $introns_dir . ${species} . "$type" . ".intron.tbl";
    fasta_2_tbl( $intron_tmp_fa, $intron_tmp_tbl );

## INTRONS LARGER THAN 0 ONLY

    my $intron_nonzero_tmp = "";
    my $my_command =
      "gawk '{if(length(\$2)>0){print \$1,\$2}}' $intron_tmp_tbl ";
    $intron_nonzero_tmp = capture($my_command);

    my $intron_nonzero_tbl =
      $work_dir . $species . "$type" . ".intron_positivelength.tbl";

    open( my $fh_FOUT, ">", "$intron_nonzero_tbl" ) or croak "Failed here";
    print $fh_FOUT "$intron_nonzero_tmp";
    close $fh_FOUT;

    print STDERR
      "intron tabular file created with introns with more than 0 nucleotides\n";
## GET LIST OF SEQUENCES WITH LENGTH >0 and EXCLUDE FROM CDS/locus_id/gff FILES SEQUENCES WITH INTRONS WITH 0 LENGTH
    my $intron_zero = "";
    $my_command =
"gawk '{if(length(\$2)==0){print \$1}}' $intron_tmp_tbl | sed 's/\\(.*\\)\\..*/\\1\\_/' | sort | uniq ";
    $intron_zero = capture($my_command);

    my $tempall_intron_zero_list =
      $work_dir . $species . "$type" . ".intron_zerolength.list";

    open( $fh_FOUT, ">", "$tempall_intron_zero_list" );
    print $fh_FOUT "$intron_zero";
    close $fh_FOUT;

    my $intron_zero2 = "";
    $my_command =
"gawk '{if(length(\$2)==0){print \$1}}' $intron_tmp_tbl | sed 's/\\(.*\\)\\..*/\\1/' | sort | uniq ";
    $intron_zero2 = capture($my_command);

    my $tempall_intron_zero_list2 =
      $work_dir . $species . "$type" . ".intron_zerolength.list2";

    open( $fh_FOUT, ">", "$tempall_intron_zero_list2" );
    print {$fh_FOUT} "$intron_zero2";
    close $fh_FOUT;

## FILTER SEQUENCES WITH 0 SIZE INTRONS FROM CDS!

    $my_command = "egrep -vf $tempall_intron_zero_list $cds_temp_tbl ";
    my $cds_nozero_tbl = capture($my_command);

    my $cds_all_nozero_tbl = $work_dir . $species . "$type" . ".cds_nozero.tbl";

    open( $fh_FOUT, ">", "$cds_all_nozero_tbl" );
    print {$fh_FOUT} "$cds_nozero_tbl";
    close $fh_FOUT;
## ENSURE LOCUSID DOES NOT CONTAIN SEQUENCES WITH 0 SIZE INTRONS
    $my_command = "egrep -vwf $tempall_intron_zero_list2 $locus_id ";
    my $locusid_nozero = capture($my_command);

    my $templocus_id_nozero =
      $work_dir . $species . "$type" . "_locus_id_nozero";
    open( $fh_FOUT, ">", "$templocus_id_nozero" );
    print {$fh_FOUT} "$locusid_nozero";
    close $fh_FOUT;
## ENSURE GFF DOES NOT CONTAIN SEQUENCES WITH 0 SIZE INTRONS

    my $gffnozero = "";
    $my_command = "egrep -vwf $tempall_intron_zero_list2 $no_dots_gff_fn ";
    $gffnozero  = capture($my_command);

    my $exCI_temp_nonzero_gff =
      $work_dir . $species . "$type" . ".non_zero.gff";

    open( $fh_FOUT, ">", "$exCI_temp_nonzero_gff" ) or croak "Failed here";
    print {$fh_FOUT} "$gffnozero";
    close $fh_FOUT;

    #    rmtree([ "$path/cds/" ]);
    #    rmtree([ "$path/intron/" ]);

## Convert sequences to protein format and check for in-frame stops
    print STDERR
"\nConvert sequences to protein format and check for in-frame stops and for proteins not starting with an M or not ending with a STOP\n\n";

## SHOWS WHERE GENETIC CODE FILE IS LOCATED AND ITS NAME

    my $tempall_protein = $work_dir . $species . "$type" . ".protein";

# $tempall_protein = translate_2_protein($genetic_code,$tempcds,$tempall_protein);
    $tempall_protein =
      translate_2_protein( $genetic_code, $cds_all_nozero_tbl,
        $tempall_protein );

    $my_command =
"gawk '{print \$2,\$1}' $tempall_protein | egrep '[A-Z]\\*[A-Z]\|^[^M]\|[^\\*] ' | gawk '{print \$2}' | wc -l";
    ## BUG reports just the number
    my $inframestops = capture($my_command);
    chomp $inframestops;

#~ #my $inframe_Xstops =`gawk '{print \$2,\$1}' $tempall_protein | egrep '[A-Z]\\*[A-Z]\|^[^M]\|[^\\*] ' | gawk '{print \$2}' | wc | gawk '{print \$1}'`;
#~ chomp $inframe_Xstops;

    print STDERR
"\n\nWe found $inframestops sequences with in-frame stop signals/not starting with a methionine or not ending with a canonical stop codon \n\n";

## IF INFRAME
    print $tempall_protein;

    if ($inframestops) {
        my $inframe = "";
        my @inframe = ();

        ## XXX may not work XXX
        $my_command =
"gawk '{print \$2,\$1}' $tempall_protein | egrep '[A-Z]\\*[A-Z]|^[^M]|[^\\*]' | gawk '{print \$2}' | sort | uniq ";
        @inframe = capture($my_command);

        foreach my $line (@inframe) {
            my (@frame) = split "_", $line;
            my $first = $frame[0];
            $inframe .= "$first\n";
        }

        my $inframe_protein =
          $work_dir . $species . "$type" . "_INframe_NoMethionine_NoSTOP";
        open( $fh_FOUT, ">", "$inframe_protein" ) or croak "Failed here";
        print {$fh_FOUT} "$inframe";
        close $fh_FOUT;
## REMOVE SEQUENCES WITH IN-FRAME STOPS FROM ORIGINAL CDS / INTRON / LOCUS_ID /GFF FILES AND PRINT NEW FILES
        print STDERR
"\nremove sequences with in-frame stop signals from cds/intron files\n\n";

        $my_command =
"sed 's/\\(.*\\)/\\1_/g' $inframe_protein | egrep -vf - $cds_all_nozero_tbl";
        my $cdstbl2 = capture($my_command);

        my $tempall_cds2 = $work_dir . $species . "$type" . ".cds_filter1.tbl";

        open( $fh_FOUT, ">", "$tempall_cds2" );
        print {$fh_FOUT} "$cdstbl2";
        close $fh_FOUT;

        my $introntbl2 = "";

        $my_command =
"sed 's/\\(.*\\)/\\1\.i/g' $inframe_protein | egrep -vf - $intron_nonzero_tbl ";
        $introntbl2 = capture($my_command);

        my $tempall_intron2 =
          $work_dir . $species . "$type" . ".intron_filter1.tbl";

        open( $fh_FOUT, ">", "$tempall_intron2" );
        print {$fh_FOUT} "$introntbl2";
        close $fh_FOUT;

        $my_command =
"sed 's/\\(.*\\)/\\1\$/g' $inframe_protein | egrep -vf - $templocus_id_nozero ";
        my $new_locus_id_filter1 = capture($my_command);

        my $templocus_id_new2 =
          $work_dir . $species . "$type" . "_locus_id_filter_noinframe";

        open( $fh_FOUT, ">", "$templocus_id_new2" );
        print {$fh_FOUT} "$new_locus_id_filter1";
        close $fh_FOUT;

        #my $gffnew = "";
        $my_command =
"sed 's/\\(.*\\)_.*/\\1\$/g' $inframe_protein | egrep -vf - $exCI_temp_nonzero_gff ";
        my $gffnew = capture($my_command);

        my $tempnewgff = $work_dir . $species . "$type" . ".noinframe.gff";

        open( $fh_FOUT, ">", "$tempnewgff" ) or croak "Failed here";
        print {$fh_FOUT} "$gffnew";
        close $fh_FOUT;

######
        return [
            $tempall_cds2, $tempall_intron2, $templocus_id_new2,
            $tempnewgff,   $inframestops
        ];

    }
    else {    ## ??? END IF THERE ARE INFRAME STOPS
        return [
            $cds_all_nozero_tbl,  $intron_nonzero_tbl,
            $templocus_id_nozero, $exCI_temp_nonzero_gff,
            0
        ];
    }    #######END ELSE IF NO SEQS  ARE INFRAME

}    #########sub extractCDSINTRON

## FUNCTION TO EXTRACT AND PROCESS SPLICE SITES AND START CODON
sub extractprocessSITES {

    my ( $no_dots_gff_fn, $locus_id ) = @_;
    my $fh_LOCID;

## SPLICE SITES
    print STDERR "\nEXTRACT START AND SPLICE SITES from transcripts\n\n";

    #print STDERR "$locus_id and $input_gff_fn\n";
    my @newsites = ();
    my $count    = 0;

    open( my $fh_LOC_sites, "<", "$locus_id" ) or croak "Failed here";
    while (<$fh_LOC_sites>) {
        my ( $genomic_id, $gene_id ) = split;

        #  print STDERR "$genomic_id,$gene_id\n";
        run("egrep -w '$gene_id\$' $no_dots_gff_fn > $tmp_dir/$gene_id.gff");

        #  print STDERR "$gene_id $input_gff_fn $tmp_dir/$gene_id.gff \n\n";
        ## POTENTIAL BUG SPLIT
        my $my_command =
"./bin/ssgff -dabeE $fastas_dir/$genomic_id $tmp_dir/$gene_id.gff > $tmp_dir/${gene_id}.all_sites";
        run($my_command);

        foreach my $site (qw(Acceptor Donor Stop Start)) {

#    print STDERR "egrep -A 1 $site $tmp_dir/${gene_id}.all_sites $sitesdir/${site}_sites.fa\n";
## POTENTIAL BUG, split command below

            run(
" egrep -A 1 $site $tmp_dir/${gene_id}.all_sites | sed -e '/--/d' -e '/^\$/d' >> $sites_dir/${site}_sites.fa"
            );
        }
        $count++;
        print STDERR "$count..";
    }    #while $fh_LOC_sites
    close $fh_LOC_sites;

    my $accesions_fa = "$sites_dir/Acceptor_sites.fa";
    my $donors_fa    = "$sites_dir/Donor_sites.fa";
    my $ATGx_fa      = "$sites_dir/Start_sites.fa";
    my $stops_fa     = "$sites_dir/Stop_sites.fa";

    my $preATGx_tbl = "$work_dir/Start_sites.tbl";
    fasta_2_tbl( $ATGx_fa, $preATGx_tbl );
    my $acceptors_tbl = "$work_dir/Acceptor_sites.tbl";
    fasta_2_tbl( $accesions_fa, $acceptors_tbl );

    #print STDERR "$acceptortbl\n";
    my $donors_tbl = "$work_dir/Donor_sites.tbl";
    fasta_2_tbl( $donors_fa, $donors_tbl );

##ADD N TO START SITES############
    ## POTENTIAL BUG
    my $ATGx_tbl = "$sites_dir" . "Start_sites_complete.tbl";
    my $my_command =
"gawk '{printf \$1\" \";for (i=1;i<=60-length(\$2);i++) printf \"n\"; print \$2}' $preATGx_tbl > $ATGx_tbl";
    run($my_command);

#`gawk '{printf \$1" ";for (i=1;i<=60-length(\$2);i++) printf "n"; print \$2}' $prestarttbl > $sites_dir/Start_sites_complete.tbl`;

#################################

    print STDERR "\n\nEliminate non-canonical donors/acceptors/starts:\n";

    #      ##EXTRACT NON CANONICAL DONORS
    #my $noncanonical     = "";
    my $generic_noncanonical = "";
    my $tot_generic_noncanonical  = "";
    my $tot_generic_canonical     = "";
    
    
    
    my $newdonortbl      = "";
    $my_command =
      "gawk '{print \$2}' $donors_tbl  | egrep -v '^[NATCGn]{31}GT' ";
    my $noncanonical = capture($my_command);

    my $donor_noncanonical_temp = $work_dir . $species . "_non_canonical_donor";
    open( my $fh_FOUT, ">", "$donor_noncanonical_temp" ) or croak "Failed here";
    print {$fh_FOUT} "$noncanonical";
    close $fh_FOUT;

    $tot_generic_noncanonical = num_of_lines_in_file($donor_noncanonical_temp);

    print STDERR
"\nThere are $tot_generic_noncanonical non-canonical donors within the training set:\n";

###########################
    if ($tot_generic_noncanonical) {    #if there are non canonical donors

        my @generic_noncanonical = ();
        open $fh_LOCID, "-|",
"egrep -wf $donor_noncanonical_temp $donors_tbl | gawk '{print \$1}' - | sort | uniq";
        while (<$fh_LOCID>) {

            push( @generic_noncanonical, "$_" );

        }
        close $fh_LOCID;

        foreach my $line (@generic_noncanonical) {

            #   my (@noncan)= split (/\.\d+:/, $line);
            my (@noncan) = split( /:/, $line );
            my $first = $noncan[0] . ":";
            $generic_noncanonical .= "$first\n";

        }

        #  unlink $donor_noncanonical_temp;

        my $noncanonical_name_tmp =
          $work_dir . $species . "_non_canonical_donor_seq_name";
        open( my $fh_FOUT, ">", "$noncanonical_name_tmp" )
          or croak "Failed here";
        print {$fh_FOUT} "$generic_noncanonical";
        close $fh_FOUT;

        open $fh_LOCID, "-|", "egrep -vf $noncanonical_name_tmp $donors_tbl";
        while (<$fh_LOCID>) {
            $newdonortbl .= $_;
        }
        close $fh_LOCID;

        my $tempcanonicaldonor = $work_dir . $species . ".canonical.donor.tbl";
        open( $fh_FOUT, ">", "$tempcanonicaldonor" ) or croak "Failed here";
        print {$fh_FOUT} "$newdonortbl";
        close $fh_FOUT;

        # unlink $noncanonical_name_tmp;

        $tot_generic_canonical = num_of_lines_in_file($tempcanonicaldonor);

        print STDERR
"\nThere are $tot_generic_canonical canonical donors within the training set:\n";

        push( @newsites, "$tempcanonicaldonor" );
        push( @newsites, "$tot_generic_noncanonical" );
    }
    else {    #if there are no non-canonical
        my $tot_generic_canonical = num_of_lines_in_file($donors_tbl);

        print STDERR
          "There are $tot_generic_canonical canonical donors within the training set:\n";
        push( @newsites, "$donors_tbl" );
        push( @newsites, "" );

    }    #if there are no non-canonical

    #      ###########################

    #      ####
    #      ##EXTRACT NON CANONICAL ACCEPTORS
    $noncanonical     = "";
    $generic_noncanonical = "";
    $tot_generic_canonical     = "";
    my $acceptor_new_tbl = "";
    #~ my $foobar_tmp     = "";

    #$my_command =
    #  "gawk '{print \$2}' $acceptortbl | egrep -v '^[NATCG]{28}AG'";

    # BUG this blows if there are no such sites...
    #$foobar_tmp = capture($my_command);

    open $fh_LOCID, "-|",
      "gawk '{print \$2}' $acceptors_tbl | egrep -v '^[NATCG]{28}AG'";
    while (<$fh_LOCID>) {
        $noncanonical .= $_;
    }
    if ( length($noncanonical) > 0 ) {
        close $fh_LOCID;
    }

    my $acceptor_noncanonical_temp =
      $work_dir . $species . "_non_canonical_acceptor";
    open( $fh_FOUT, ">", "$acceptor_noncanonical_temp" ) or croak "Failed here";
    print {$fh_FOUT} "$noncanonical";
    close $fh_FOUT;

    $tot_generic_noncanonical = num_of_lines_in_file($acceptor_noncanonical_temp);

    print STDERR
"\nThere are $tot_generic_noncanonical non-canonical acceptors within the training set:\n";
###########################
    if ($tot_generic_noncanonical) {    #if there are non-canonical acceptors

        my @generic_noncanonical = ();
        open $fh_LOCID, "-|",
"egrep -f $acceptor_noncanonical_temp $acceptors_tbl | gawk '{print \$1}' - | sort | uniq ";
        while (<$fh_LOCID>) {
            push( @generic_noncanonical, "$_" );
        }

        close $fh_LOCID;

        foreach my $line (@generic_noncanonical) {
            my (@noncan) = split( /:/, $line );
            my $first = $noncan[0] . ":";
            $generic_noncanonical .= "$first\n";

        }

        #~ unlink $acceptor_noncanonical_temp;

        my $noncanonical_name_tmp =
          $work_dir . $species . "_non_canonical_acceptor_seq_name";

        open( my $fh_FOUT, ">", "$noncanonical_name_tmp" );
        print {$fh_FOUT} "$generic_noncanonical";
        close $fh_FOUT;

        open $fh_LOCID, "-|", "egrep -vf $noncanonical_name_tmp $acceptors_tbl";
        while (<$fh_LOCID>) {
            $acceptor_new_tbl .= $_;
        }
        close $fh_LOCID;

        #unlink $noncanonical_name_tmp;

        my $acceptor_canonical_temp =
          $work_dir . $species . ".canonical.acceptor.tbl";
        open( $fh_FOUT, ">", "$acceptor_canonical_temp" ) or croak "Failed here";
        print {$fh_FOUT} "$acceptor_new_tbl";
        close $fh_FOUT;

        #unlink $noncanonical_name_tmp;

        my $tot_generic_canonical = num_of_lines_in_file($acceptor_canonical_temp);

        print STDERR
"\nThere are $tot_generic_canonical canonical acceptors within the training set:\n";

        push( @newsites, "$acceptor_canonical_temp" );
        push( @newsites, "$tot_generic_noncanonical" );

    }
    else {    #if there are only canonical use initial file list
        my $tot_generic_canonical = num_of_lines_in_file($acceptors_tbl);

        print STDERR
"There are $tot_generic_canonical canonical acceptors within the training set:\n";
        push( @newsites, "$acceptors_tbl" );
        push( @newsites, "0" );
    }    #if there are only canonical use initial file list

    #      ###########################

    #      ###
    #      ##EXTRACT NON CANONICAL STARTS

    $noncanonical     = "";
    $generic_noncanonical = "";
    $tot_generic_canonical     = "";
    my $newstarttbl = "";

    open $fh_LOCID, "-|",
      "gawk '{print \$2}' $ATGx_tbl | egrep -v '^[NATCG]{30}ATG' ";
    while (<$fh_LOCID>) {
        $noncanonical .= $_;
    }
    if ( length($noncanonical) > 0 ) {
        close $fh_LOCID;
    }

    #close $fh_LOCID;

    my $tempstartnoncanonical = $work_dir . $species . "_non_canonical_start";
    open( $fh_FOUT, ">", "$tempstartnoncanonical" ) or croak "Failed here";
    print {$fh_FOUT} "$noncanonical";
    close $fh_FOUT;

    $tot_generic_noncanonical = num_of_lines_in_file($tempstartnoncanonical);

    print STDERR
"\nThere are $tot_generic_noncanonical non-canonical starts within the training set:\n";
###########################

    if ($tot_generic_noncanonical) {    #if there are non-canonical starts

        my @generic_noncanonical = ();
        open $fh_LOCID, "-|",
"egrep -wf $tempstartnoncanonical $ATGx_tbl | gawk '{print \$1}' - | sort | uniq ";
        while (<$fh_LOCID>) {
            push( @generic_noncanonical, "$_" );
        }
        close $fh_LOCID;

        foreach my $line (@generic_noncanonical) {
            my (@noncan) = split( /:/, $line );
            my $first = $noncan[0] . ":";
            $generic_noncanonical .= "$first\n";

        }

        unlink $tempstartnoncanonical;

        my $noncanonical_name_tmp =
          $work_dir . $species . "_non_canonical_start_seq_name";
        open( $fh_FOUT, ">", "$noncanonical_name_tmp" );
        print {$fh_FOUT} "$generic_noncanonical";
        close $fh_FOUT;

        open( $fh_LOCID, "-|", "egrep -vf $noncanonical_name_tmp $ATGx_tbl " );
        while (<$fh_LOCID>) {
            $newstarttbl .= $_;
        }
        close $fh_LOCID;

        # unlink $noncanonical_name_tmp;

        my $tempcanonicalstart = $work_dir . $species . ".canonical.start.tbl";
        open( $fh_FOUT, ">", "$tempcanonicalstart" ) or croak "Failed here";
        print {$fh_FOUT} "$newstarttbl";
        close $fh_FOUT;

        #unlink $noncanonical_name_tmp;
        my $tot_generic_canonical = num_of_lines_in_file($tempcanonicalstart);

        print STDERR
"\nThere are $tot_generic_canonical canonical starts within the training set:\n";

        push( @newsites, "$tempcanonicalstart" );
        push( @newsites, "$tot_generic_noncanonical" );

    }
    else {
        my $tot_generic_canonical = num_of_lines_in_file($ATGx_tbl);

        print STDERR
"\nThere are $tot_generic_canonical canonical starts within the training set:\n";
        push( @newsites, "$ATGx_tbl" );
        push( @newsites, "0" );

    }    #if there are only canonical starts
###########################

    return \@newsites;

}    #subectractprocesssites

## FUNCTION TO OBTAIN MARKOV MODELS CORRESPONDING TO THE CODING POTENTIAL
sub deriveCodingPotential {
    my ( $cds, $intron ) = @_;

    my $markov_mod_A = "";
    my $markov_mod_B = "";

    my $my_command = "gawk '{ l=length(\$2); L+=l;} END{ print L;}' $cds";
    my $total_codingbases = capture($my_command);
    chomp $total_codingbases;

    $my_command = "gawk '{ l=length(\$2); L+=l;} END{ print L;}' $intron ";
    my $total_noncodingbases = capture($my_command);
    chomp $total_noncodingbases;

    print STDERR
"There are $total_codingbases coding bases and $total_noncodingbases non-coding bases on this training set:\n";

    if (
        (
               $total_codingbases > $coding_bp_limit_A
            && $total_noncodingbases > $non_coding_bp_limit_A
        )
        || (   $total_codingbases > $coding_bp_limit_B
            && $total_noncodingbases > $non_coding_bp_limit_B )
        || (   $total_noncodingbases > 35_000
            && $total_codingbases > ( 25 * $total_noncodingbases ) )
      )
    {
        $markov_mod_A = 5;
        $markov_mod_B = 4;
        print STDERR
          "Deriving a markov model of order $markov_mod_A OPTION_1\n";

    }
    else {
        $markov_mod_A = 5;
        $markov_mod_B = 4;
        print STDERR
          "Deriving a markov model of order $markov_mod_A  OPTION_2\n";
    }

    open( my $fh_INTRONS, "<", "$intron" ) or croak "Failed here";
    my @introndois = ();
    while (<$fh_INTRONS>) {
        my @i = split;

        #print STDERR "SECOND FIELD $i[2]";
        push @introndois, $i[1];
    }
    close $fh_INTRONS;

    open( my $fh_CDSes, "<", "$cds" ) or croak "Failed here";
    my @coding = ();
    while (<$fh_CDSes>) {
        my @c = split;
        push @coding, $c[1];
    }
    close $fh_CDSes;

    print STDERR "Intron model\n markov: ($markov_mod_B)";

    my $intron_initial =
      geneidCEGMA::SequenceModel->new( 'intron', 'FREQ', $markov_mod_B,
        \@introndois, 10, 0 );

    my $intron_transition =
      geneidCEGMA::SequenceModel->new( 'intron', 'MM', $markov_mod_A,
        \@introndois, 10, 0 );

    print STDERR "Coding model\n";

    my $coding_initial =
      geneidCEGMA::SequenceModel->new( 'coding', 'FREQ', $markov_mod_A - 1,
        \@coding, 0.25, 2 );

    my $coding_transition =
      geneidCEGMA::SequenceModel->new( 'coding', 'MM', $markov_mod_A, \@coding,
        0.25, 2 );

    my $initial_logs =
      geneidCEGMA::log_ratio( $coding_initial, $intron_initial );

    my $transition_logs =
      geneidCEGMA::log_ratio( $coding_transition, $intron_transition );

    geneidCEGMA::write_log( $initial_logs, "$work_dir/coding.initial.5.logs" );

    geneidCEGMA::write_log( $transition_logs,
        "$work_dir/coding.transition.5.logs" );

    open( my $fh_PROFILE_1, "<", "$work_dir/coding.initial.5.logs" )
      or croak "Failed here";
    my @profileinit = ();
    while (<$fh_PROFILE_1>) {
        last if m/^\s/;
        last if m/^[^ACGTacgt]/;
        next if m/^#/;
        chomp;
        my @g = split;
        push @profileinit, \@g;
    }
    close $fh_PROFILE_1;

    open( my $fh_PROFILE_2, "<", "$work_dir/coding.transition.5.logs" )
      or croak "Failed here";
    my @profiletran = ();
    while (<$fh_PROFILE_2>) {
        last if m/^\s/;
        last if m/^[^ACGTacgt]/;
        next if m/^#/;
        chomp;
        my @g = split;
        push @profiletran, \@g;
    }
    close $fh_PROFILE_2;

    return [
        \@profileinit,         \@profiletran, $total_codingbases,
        $total_noncodingbases, $markov_mod_A
    ];

}    #derive coding potential

## PROCESS SEQUENCES FUNCTION ( FLANKED GENE MODELS OBTAINED FOR OPTIMIZATION)
sub processSequences4Optimization {

    my ( $no_dots_gff_fn, $type, $run_contig_opt_flag ) = @_;

    my $pso_out_tbl = "";
    my $pso_gp_tbl  = "";
    my $gp_from_gff = "";

    #~ my $gp_fasta    = ""; #unused
    my $pso_gp_gff = "";
    my $my_command = "";

    #my $work_dir;

    open( my $fh_LOCID, "-|", "./bin/gff2gp.awk $no_dots_gff_fn | sort -k 1 " );
    while (<$fh_LOCID>) {

        $gp_from_gff .= $_;
    }
    close $fh_LOCID;

    my $pso_tmp_gp_from_gff = $work_dir . $species . $type . ".gp";

    open( my $fh_FOUT, ">", "$pso_tmp_gp_from_gff" );
    print {$fh_FOUT} "$gp_from_gff";
    close $fh_FOUT;
    print STDERR
"BEFORE GETGENES: $fastas_dir, $pso_tmp_gp_from_gff, $work_dir/, $pso_out_tbl\n";

    my $gp_Xgetgenes_tmp_pre_tbl =
      GetGenes( $fastas_dir, $pso_tmp_gp_from_gff, $work_dir, $pso_out_tbl );
    print STDERR "PRETBL AFTER GETGENES: $gp_Xgetgenes_tmp_pre_tbl \n";

    print STDERR
"\nGet sequences of 400-nt flanked sequences in tabular and gff formats\n";

#~ my $seq4Optimization_temp_1_fn =  "$work_dir/processSequences4Optimization_temp1.txt";
#~ $my_command =  "gawk 'BEGIN{b=\"x\"}{if (\$1!=b){printf \"\\n\%s \",\$1}printf \"\%s\",\$6;b=\$1}END{printf \"\\n\"}' $gp_Xgetgenes_tmp_pre_tbl > $seq4Optimization_temp_1_fn";
#~ run($my_command);
#~ $my_command =  "sort seq4Optimization_temp_1_fn | sed '/^\$/d' - | gawk '{print \$1, toupper(\$2)}' - |";
#open( $fh_LOCID, $my_command );
    open( $fh_LOCID, "-|",
"gawk 'BEGIN{b=\"x\"}{if (\$1!=b){printf \"\\n\%s \",\$1}printf \"\%s\",\$6;b=\$1}END{printf \"\\n\"}' $gp_Xgetgenes_tmp_pre_tbl | sort | sed '/^\$/d' - | gawk '{print \$1, toupper(\$2)}' - "
    );

    while (<$fh_LOCID>) {
        $pso_gp_tbl .= $_;
    }
    close $fh_LOCID;

    my $tempgp_tbl = $work_dir . $species . $type . ".gp.tbl";
    open( $fh_FOUT, ">", "$tempgp_tbl" ) or croak "Failed here";
    print {$fh_FOUT} "$pso_gp_tbl";
    close $fh_FOUT;

    open( $fh_LOCID, "-|",
"gawk 'BEGIN{OFS=\"\\t\";pos=1;b=\"x\"}{if (\$1!=b){pos=1}; print \$1,\"annotations\",\$3,pos,pos+\$5-1,\"\.\",\"+\",\"\.\",\$1\$2; pos+=\$5;b=\$1 }' $gp_Xgetgenes_tmp_pre_tbl | egrep -v '(Intron|Utr)' - "
    );
    while (<$fh_LOCID>) {
        $pso_gp_gff .= $_;
    }
    close $fh_LOCID;

    my $tempgp_gff = $work_dir . $species . $type . ".gp.gff";
    open( $fh_FOUT, ">", "$tempgp_gff" ) or croak "Failed here";
    print {$fh_FOUT} "$pso_gp_gff";
    close $fh_FOUT;

    print STDERR "DONE\n";

    print STDERR
      "\nGet sequences of 400-nt flanked sequences in multi-fasta format\n";

    my $tempgp_fa = $work_dir . $species . $type . ".gp.fa";

    $tempgp_fa = TblToFasta( $tempgp_tbl, $tempgp_fa );

    ## BUG keep intermediates
    #~ unlink $gp_Xgetgenes_tmp_pre_tbl;

    print STDERR "\nSet up files for optimization\n\n";
    $my_command = "gawk '{print \$1,length(\$2)}' $tempgp_tbl | sort -k1,1 ";
    my $seqslenggp = capture($my_command);

    #  ` gawk '{print \$1,length(\$2)}' $tempgp_tbl | sort -k1,1 `;    ##XX

    my $tempseqlen = $work_dir . $species . $type . ".gp_cds_length";
    open( $fh_FOUT, ">", "$tempseqlen" ) or croak "Failed here";
    print {$fh_FOUT} "$seqslenggp";
    close $fh_FOUT;

    my $cdsgp = "";
    open( $fh_LOCID, "-|",
"./bin/gff2cds.awk source=\"annotations\" $tempgp_gff | sort -k1,1 | join $tempseqlen - "
    );
    while (<$fh_LOCID>) {
        $cdsgp .= $_;
    }
    close $fh_LOCID;

    my $tempcdsgp = $work_dir . $species . $type . ".cds_gp";
    open( $fh_FOUT, ">", "$tempcdsgp" ) or croak "Failed here";
    print {$fh_FOUT} "$cdsgp";
    close $fh_FOUT;

    my $pso_gp_eval_gff = "";
    open( $fh_LOCID, "-|",
"gawk 'BEGIN{while (getline<ARGV[1]>0){len[\$1]=\$2;};ARGV[1]=\"\";OFS=\"\\t\";}{if (NR==1) {ant=\$1;print \$1,\$2,\"Sequence\",1,len[ant],\"\.\",\"\.\",\"\.\",\"\.\"};if (\$1!=ant) {print \"\#\$\";ant=\$1;print \$1,\$2,\"Sequence\",1,len[ant],\"\.\",\"\.\",\"\.\",\"\.\"}; print }' $tempcdsgp $tempgp_gff "
    );
    while (<$fh_LOCID>) {
        $pso_gp_eval_gff .= $_;

    }
    close $fh_LOCID;

    my $tempevalgpgff = $work_dir . $species . $type . ".gp_eval_gff";
    open( $fh_FOUT, ">", "$tempevalgpgff" ) or croak "Failed here";
    print {$fh_FOUT} "$pso_gp_eval_gff";
    close $fh_FOUT;

    if ($run_contig_opt_flag) {

        #my $tempgp_fa = $species.$type.".gp.fa";

        #$tempgp_fa = FastaToTbl($tempgp_tbl,$tempgp_fa);

        my @gp_tabular = split( /\n/, $pso_gp_tbl );
        my $seq = "";
        foreach my $line (@gp_tabular) {
            chomp $line;
            my @f = split " ", $line;
            $seq .= $f[1];
        }
        my $lengp       = length($seq);
        my $foldedseqgp = fold4fasta($seq);
        my $tempfastagpcontig =
          $work_dir . $species . $type . ".combined.gp.fa";
        open( $fh_FOUT, ">", "$tempfastagpcontig" ) or croak "Failed here";
        print {$fh_FOUT} ">$species\n$foldedseqgp\n";
        close $fh_FOUT;

        my $temptabulargpcontig =
          $work_dir . $species . $type . ".combined.gp.tbl";
        open( $fh_FOUT, ">", "$temptabulargpcontig" ) or croak "Failed here";
        print {$fh_FOUT} "$species\t$seq\n";
        close $fh_FOUT;
        $my_command = "gawk '{print \$1,length(\$2)}' $temptabulargpcontig";
        my $seqslengcontiggp = capture($my_command);

        #` gawk '{print \$1,length(\$2)}' $temptabulargpcontig `;

        my $tempseqlencontig =
          $work_dir . $species . $type . ".gp_cds_contig_length";
        open( $fh_FOUT, ">", "$tempseqlencontig" ) or croak "Failed here";
        print {$fh_FOUT} "$seqslengcontiggp";
        close $fh_FOUT;

        my $gpcontig = "";
        open( $fh_LOCID, "-|",
"./bin/multiple_annot2one.awk species=$species leng=$lengp $tempcdsgp "
        );
        while (<$fh_LOCID>) {

            $gpcontig .= $_;
        }
        close $fh_LOCID;

        my $tempgff2gpcontig = $work_dir . $species . $type . ".contig.gp.cds";

        open( $fh_FOUT, ">", "$tempgff2gpcontig" ) or croak "Failed here";
        print {$fh_FOUT} "$gpcontig";
        close $fh_FOUT;

        my $cds2gffcontig = "";
        open( $fh_LOCID, "-|",
"./bin/cds2gff.awk $tempgff2gpcontig | gawk 'BEGIN{OFS=\"\\t\";}{if (NR==1){print \"$species\",\"annotations\",\"Sequence\",\"1\",$lengp,\".\",\".\",\".\",\".\";print}else {print}}' - "
        );
        while (<$fh_LOCID>) {
            $cds2gffcontig .= $_;
        }
        close $fh_LOCID;

        my $tempgp_cdsgff_contig_eval =
          $work_dir . $species . $type . ".cds_gp_contig.eval.gff";
        open( $fh_FOUT, ">", "$tempgp_cdsgff_contig_eval" )
          or croak "Failed here";
        print {$fh_FOUT} "$cds2gffcontig";
        close $fh_FOUT;

        return [
            $tempgp_cdsgff_contig_eval, $tempfastagpcontig,
            $temptabulargpcontig,       $tempseqlencontig
        ];

    }
    elsif ( !$run_contig_opt_flag ) {
        print STDERR "L1803, NOT CONTIG OPT\n";
        return [ $tempevalgpgff, $tempgp_fa, $tempgp_tbl, $tempseqlen ];
    }

}    #processSequences optimization

## GETGENES FUNCTION: EXTRACT FLANKED SEQUENCES FROM GENE MODELS FOR LATER OPTIMIZATION
sub GetGenes {

    my ( $path2gpath, $genes_fn_X, $work_dir, $gp_out_tbl ) = @_;

#print STDERR "IN FUNCTION: $path2gpath : $genes_fn_X : $path : OUT: $gp_out_tbl\n\n";

    my $nonred     = 0;
    my $onlynonred = 0;
    my $prevgene   = "x";
    my $prevchro   = "x";
    my $trail      = "";

    #my %genenames; unused var
    $gp_out_tbl = "$work_dir/gp_out_X.tbl";

    #~ chomp($path2gpath);
    #~ chomp($genes_fn_X);
    #~ chomp($work_dir);
    #~ chomp($gp_out_tbl);

    open( my $fh_REFGENE,   "<", "$genes_fn_X" ) or croak "Failed here";
    open( my $fh_OUT_tblgb, ">", "$gp_out_tbl" ) or croak "Failed here";
    while (<$fh_REFGENE>) {

m/([\w\-\.:]+)\s+([\w\.\-:]+)\s+([\+\-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([^\n]+)/;

        my $name      = $1;
        my $chro      = $2;
        my $stra      = $3;
        my $txSt      = $4;
        my $txEn      = $5;
        my $cdsS      = $6;
        my $cdsE      = $7;
        my $exoC      = $8;
        my @exon      = ( $9 =~ m/(\d+)/g );
        my $cdsLe     = $cdsE - $cdsS;
        my $txLe      = $txEn - $txSt;
        my $cdsoffset = $cdsS - $txSt;
        my $redundant = 0;
        my $i         = 0;                     # exon counter
        my $j         = 0;
        my $call      = "";
        my $sub_seq   = "";
        my $genomic   = "";

        #my @tabular = ();

        if ( !$onlynonred || ( $onlynonred && !$redundant ) ) {
            open( my $fh_FLEN, "<", "$path2gpath${chro}_len" )
              or croak "Failed here";
            my $line = <$fh_FLEN>;
            chomp $line;
            my @le = split " ", $line;
            close $fh_FLEN;
###added code
            my $chrotmptbl = "$tmp_dir/tmp.tbl";
            $chrotmptbl = fasta_2_tbl( $path2gpath . $chro, $chrotmptbl );

            #print STDERR "FATOTBL: $path2gpath"."$chro\n";
            open( my $fh_IN, "<", "$chrotmptbl" ) or croak "Failed here";
            my @tabular = ();
            while (<$fh_IN>) {

                #chomp;
                #print STDERR "$_";
                push @tabular, "$_";
            }
            close $fh_IN;

            #  print STDERR "\nGP: @tabular\n";
            # my @tabular = ` FastaToTbl $path2gpath$chro `;
            my $sub_seq = "";

            #my $sublen = 0;
            foreach my $line (@tabular) {
                chomp $line;
                my @f = split " ", $line;

                #print STDERR "$f[0]\n";
                $sub_seq .= $f[1];

            }
## DEBUG added

            if ( $le[1] < $txEn ) {

                my $newlen = $le[1];
                $genomic = substr( $sub_seq, $txSt, ( $newlen - $txSt ) );

            }
            elsif ( $le[1] >= $txEn ) {
                $genomic = substr( $sub_seq, $txSt, $txLe );

            }

            # my $genomic = `$call`;
            # my $genomicLe = length($genomic);
            my $genomicLe = length($genomic);
            my $cds_seq   = "";

            if ( $genomicLe == 0 ) {
                print STDERR "getgenes: gene of 0 length ($name), $call\n";
                next;
            }

            #  if ($genomicLe != $newlen) {
            #    print STDERR "getgenes: length mismatch ($name)\n";
            #    next;
            #  }

            for ( $i = 0 ; $i < $exoC ; $i++ ) {

                my $utrB = 0;
                my $utrA = 0;
                my $utrS = 0;
                my $utrL = 0;
                my $exSt = $exon[$i] - $cdsS;
                my $exLe = $exon[ $i + $exoC ] - $exon[$i];
                my $exTy = "Internal";

                if ( $exSt + $exLe > 0 && $exSt < $cdsLe ) {    # cds

                    if ( $exSt <= 0 || $i == 0 ) {
                        if ( $stra eq '+' ) {
                            $exTy = "First";
                        }
                        else {
                            $exTy = "Terminal";
                        }
                    }

                    if ( $exSt + $exLe >= $cdsLe || $i == $exoC - 1 ) {
                        if ( $stra eq '+' ) {
                            $exTy = "Terminal";
                        }
                        else {
                            $exTy = "First";
                        }
                    }

                    if ( $exSt <= 0 && $exSt + $exLe >= $cdsLe ) {
                        $exTy = "Single";
                    }

                    if ( $exSt < 0 ) {
                        $utrB = 1;
                        $utrS = $exSt;
                        $utrL = abs($exSt);
                        $exLe = $exLe - abs($exSt);
                        $exSt = 0;
                    }

                    if ( $exSt + $exLe > $cdsLe ) {
                        $utrA = 1;
                        $utrS = $cdsLe;
                        $utrL = $exLe - ( $cdsLe - $exSt );
                        $exLe = $cdsLe - $exSt;
                    }

                    my $iex;
                    my $seq = substr( $genomic, $exSt + $cdsoffset, $exLe );

                    $seq = lc($seq);

                    if ( $stra eq '+' ) {    # forward

                        if ($utrB) {
                            my $iutr = $i + 1;
                            my $utrs =
                              substr( $genomic, $utrS + $cdsoffset, $utrL );

                            $utrs    = lc($utrs);
                            $cds_seq = $cds_seq
                              . "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n";
                        }

                        $iex     = $i + 1;
                        $cds_seq = $cds_seq
                          . "$name\t$chro\t$exTy\t$iex\t$exLe\t$seq\t$exon[$i]\t$exon[$i+$exoC]\n";

                        if ($utrA) {
                            my $iutr = $i + 1;
                            my $utrs =
                              substr( $genomic, $utrS + $cdsoffset, $utrL );

                            $utrs    = lc($utrs);
                            $cds_seq = $cds_seq
                              . "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n";
                        }

                    }
                    else {    # reverse

                        if ($utrB) {
                            my $iutr = $exoC - $i;
                            my $utrs =
                              substr( $genomic, $utrS + $cdsoffset, $utrL );

                            $utrs = lc($utrs);
                            $utrs =~ tr/acgt/tgca/;
                            $utrs = reverse($utrs);
                            $cds_seq =
                              "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n"
                              . $cds_seq;
                        }

                        $iex = $exoC - $i;
                        $seq =~ tr/acgt/tgca/;
                        $seq = reverse($seq);
                        $cds_seq =
"$name\t$chro\t$exTy\t$iex\t$exLe\t$seq\t$exon[$i+$exoC]\t$exon[$i]\n"
                          . $cds_seq;

                        if ($utrA) {
                            my $iutr = $exoC - $i;
                            my $utrs =
                              substr( $genomic, $utrS + $cdsoffset, $utrL );

                            $utrs = lc($utrs);
                            $utrs =~ tr/acgt/tgca/;
                            $utrs = reverse($utrs);
                            $cds_seq =
                              "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n"
                              . $cds_seq;
                        }

                    }

                    if (
                        $exTy ne "Single"
                        && (   ( $exTy ne "Terminal" && $stra eq '+' )
                            || ( $exTy ne "First" && $stra eq '-' ) )
                      )
                    {

                        my $inSt = $exon[ $i + $exoC ] - $cdsS;
                        my $inLe = $exon[ $i + 1 ] - $exon[ $i + $exoC ];

                        if ( $inSt + $inLe > 0 && $inSt < $cdsLe ) {

                            if ( $inSt < 0 ) {
                                print "getgenes: intron out of range! (1)\n";
                                exit(1);
                            }

                            if ( $inSt + $inLe > $cdsLe ) {
                                print "getgenes: intron out of range! (2)\n";
                                exit(1);
                            }

                            $seq =
                              substr( $genomic, $inSt + $cdsoffset, $inLe );
                            $seq = "\L$seq";

                            my $iIn;

                            if ( $stra eq '+' ) {    # forward
                                $iIn     = $j + 1;
                                $cds_seq = $cds_seq
                                  . "$name\t$chro\tIntron\t$iIn\t$inLe\t$seq\n";
                            }
                            else {
                                $iIn = $exoC - $j - 1;
                                $seq =~ tr/acgt/tgca/;
                                $seq = reverse($seq);
                                $cds_seq =
                                  "$name\t$chro\tIntron\t$iIn\t$inLe\t$seq\n"
                                  . $cds_seq;
                            }

                        }
                        else {
                            print "getgenes.pl: intron out of range! (3)\n";
                            if ( $inSt + $inLe <= 0 ) {
                                print "getgenes.pl: intron in 5' UTR\n";
                            }
                            else {
                                print "getgenes.pl: intron in 3' UTR\n";
                            }

                            exit(1);
                        }

                        $j = $j + 1;
                    }

                }
                else {    # UTRs

                    $exSt = $exon[$i] - $txSt;
                    $exLe = $exon[ $i + $exoC ] - $exon[$i];

                    my $utrs = substr( $genomic, $exSt, $exLe );

                    if ( $stra eq '+' ) {    # forward
                        my $iutr = $i + 1;

                        $utrs = lc($utrs);
                        $cds_seq =
                          $cds_seq . "$name\t$chro\tUtr\t$iutr\t$exLe\t$utrs\n";
                    }
                    else {                   # reverse
                        my $iutr = $exoC - $i;

                        $utrs = lc($utrs);
                        $utrs =~ tr/acgt/tgca/;
                        $utrs = reverse($utrs);
                        $cds_seq =
                          "$name\t$chro\tUtr\t$iutr\t$exLe\t$utrs\n" . $cds_seq;
                    }

                }
            }

            print {$fh_OUT_tblgb} $cds_seq;

        }
        elsif ($onlynonred) {
            print STDERR "$name\n";
        }

    }
    close $fh_OUT_tblgb;
    close $fh_REFGENE;

    return $gp_out_tbl;

}    #getgenes 

## GET BACKGROUND SEQUENCES (ex. 62 Kmer) used to compute log likelihoods

sub BitScoreGraph {

    my ( $infooutput, $info_thresh, $offset ) = @_;
    print STDERR "bitscoregraph input:  $infooutput, $info_thresh, $offset";
    my @info = ( $offset - 1, $offset + 1 );
    my @fields;

    open( my $fh_INFO, "<", "$infooutput" ) or croak "Failed here";
    while (<$fh_INFO>) {
        next if m/^#/;
        last if m/^\s/;
        last if m/^[^\d]/;
        chomp;
        print STDERR "QQQ prefields: $_";
        
        @fields = split;
        printf STDERR "%2s %2.2f %s",
          ##
          ( $fields[0], $fields[1], "=" x int( $fields[1] * 30 ) );
        if ( $fields[1] > $info_thresh ) {
            push( @info, $fields[0] );
        }
        print STDERR "\n";
    }
    close $fh_INFO;
    print STDERR "\n BitScoreGraph \n";

    my @sortedinfo = sort numerically @info;
    my $start      = ( shift @sortedinfo );
    if ( $start < 1 ) {
        $start = 1;
    }
    my $end = pop @sortedinfo;

    return ( $start, $end );
}    #end BitScoreGraph

## GETKMATRIX FUNCTION (Splice site an Start codon PWMs)
sub get_K_matrix {
    ## BUG do not relay on  jacknife/branch etc. use some: type ??
    #was my our?
    my ( $true_kmers_tbl, $backgrnd_kmers_tbl, $order, $offset,
        $donor, $accept, $ATG, $branch, $start, $end, $jacknife )
      = @_;

    my $matrix_type = "";
    if ($donor) {
        $matrix_type = 'donor';
    }
    if ($accept) {
        $matrix_type = 'acceptor';
    }

    if ($ATG) {
        $matrix_type = 'ATG';
    }

    #~ if ($branch) {
    #~ $matrix_type = 'branch';
    #~ }
    my $original_offset = $offset;
    my @prof            = ();
    my $tempinfolog;
    my @orders  = (qw(order-0 di tri order-4 order-5 order-6 order-7 order-8));
    my $ordname = $orders[$order];
    my $sort    = "sort -n";
    if ( $order > 1 ) {
        $sort = "sort -k1,1n -k2,2";
    }

    #    my @info = ($offset-1,$offset+1);
    my $prof_len    = 0;
    my $info_thresh = "";    #bits

    ## BUG?
    my $true_seq_name = $true_kmers_tbl;
    $true_seq_name =~ s/\.tbl$//;
    my $false_seq_name = $backgrnd_kmers_tbl;
    $false_seq_name =~ s/\.tbl$//;
    print "\n XXX L2182 true_seq_name $false_seq_name \n";

## Open true sequences
    #    print STDERR "$true_kmers_tbl (true)\n";
    open( my $fh_true_seq, "<", "$true_kmers_tbl" ) or croak "Failed here";
    $_ = <$fh_true_seq>;
    my @t   = split;
    my $len = length( $t[1] );
    close $fh_true_seq;

## Open false (background???) sequences
    #    print STDERR "$backgrnd_kmers_tbl (false)\n";
    open( my $fh_FALSE_SEQ, "<", "$backgrnd_kmers_tbl" )
      or croak "Couldn't open $backgrnd_kmers_tbl: $OS_ERROR \n";
    $_ = <$fh_FALSE_SEQ>;
    my @f    = split;
    my $len2 = length( $f[1] );
    close $fh_FALSE_SEQ;

    #    die "$len != $len2\n" if $len != $len2;
    my $true_seq_freq_fn  = $work_dir . basename($true_seq_name) . ".freq";
    my $false_seq_freq_fn = $work_dir . basename($false_seq_name) . ".freq";

#my $subtracted_true_false_freq_fn = $work_dir . basename($true_seq_name) . "_" . basename($false_seq_name).freq_subtr";
    my $my_freq_subtract_fn =
        $work_dir
      . basename($true_seq_name) . "_"
      . basename($false_seq_name)
      . ".information";

    run("./bin/frequency.py 1 $true_kmers_tbl  >  $true_seq_freq_fn");
    run("./bin/frequency.py 1 $backgrnd_kmers_tbl >  $false_seq_freq_fn");

    my $my_command_A =
      "./bin/information.py  $true_seq_freq_fn $false_seq_freq_fn ";

#~ my $my_command_B =
#~ "| gawk 'NF==2 && \$1<=$my_freq_field_limit_1 && \$1>=$my_freq_field_limit_2'";

    my $my_command_C = " > $my_freq_subtract_fn ";

    #my $my_command = $my_command_A . $my_command_B;
    my $my_command = $my_command_A . $my_command_C;

    say "\n $my_command \n";
    run($my_command);
    $tempinfolog = $my_freq_subtract_fn;

    ## True_False req_matrix_fn
    my $my_True_freq_matrix_fn =
      $work_dir . basename($true_seq_name) . "_" . "$ordname.matrix";
    ## False_True req_matrix_fn
    my $my_False_freq_matrix_fn =
      $work_dir . basename($false_seq_name) . "_" . "$ordname.matrix";

    ## True logratio req_matrix_fn
    my $my_T_generic_logratio_freq_matrix_fn =
      $work_dir . basename($true_seq_name) . "$matrix_type.log.$ordname.matrix";

    if ( !$order ) {
        $my_command =
"gawk -f ./bin/logratio_zero_order.awk $false_seq_freq_fn $true_seq_freq_fn > $my_T_generic_logratio_freq_matrix_fn";
        say "\n $my_command \n";
        run($my_command);

    }
    else {
        $my_command =
"gawk -f ./bin/Getkmatrix.awk $order $len $true_kmers_tbl | $sort > $my_True_freq_matrix_fn";
        say "\n $my_command \n";
        run($my_command);

#~ run(
#~ " gawk -f ./bin/Getkmatrix.awk $order $len $true_kmers_tbl | $sort > $true_seq_name.$ordname-matrix"
#~ );

        $my_command =
"gawk -f ./bin/Getkmatrix.awk $order $len2 $backgrnd_kmers_tbl | $sort > $my_False_freq_matrix_fn ";
        say "\n $my_command \n";
        run($my_command);

        $my_command =
"gawk -f ./bin/logratio_kmatrix.awk $my_False_freq_matrix_fn $my_True_freq_matrix_fn > $my_T_generic_logratio_freq_matrix_fn ";
        say "\n $my_command \n";
        run($my_command);

    }

    #need to check output and then go on
## draw bit score bar graph function (nested, local)

    #~ local *BitScoreGraph = sub {

    #~ my ( $infooutput, $info_thresh, $offset ) = @_;
    #~ my @info = ( $offset - 1, $offset + 1 );
    #~ my @fields;
    #~ open( my $fh_INFO, "<", "$infooutput" ) or croak "Failed here";
    #~ while (<$fh_INFO>) {
    #~ next if m/^#/;
    #~ last if m/^\s/;
    #~ last if m/^[^\d]/;
    #~ chomp;
    #~ @fields = split;
    #~ printf STDERR "%2s %2.2f %s",
    #~ ##
    #~ ( $fields[0], $fields[1], "=" x int( $fields[1] * 30 ) );
    #~ if ( $fields[1] > $info_thresh ) {
    #~ push( @info, $fields[0] );
    #~ }
    #~ print STDERR "\n";
    #~ }
    #~ close $fh_INFO;
    #~ print STDERR "\n BitScoreGraph \n";

    #~ my @sortedinfo = sort numerically @info;
    #~ my $start      = ( shift @sortedinfo );
    #~ if ( $start < 1 ) {
    #~ $start = 1;
    #~ }
    #~ my $end = pop @sortedinfo;

    #~ return ( $start, $end );
    #~ };    #end BitScoreGraph

    ## TODO simplify

    $jacknife = 0;

    #~ if ( !$jacknife ) {
    #~ print STDERR "Information content profile\n";
    #~ }

    my %my_info_thresholds = (
        donor    => 0.15,
        acceptor => 0.04,
        ATG      => 0.15,
        branch   => 0.30,
    );
    print STDERR "\n L2917: $matrix_type, $my_info_thresholds{$matrix_type}\n";

    my $my_info_thresh = $my_info_thresholds{$matrix_type};

    #say "\n matrix sub: $matrix_type, $my_info_thresh \n";
    ( $start, $end ) = BitScoreGraph( $tempinfolog, $my_info_thresh, $offset );
    print STDERR "\n L2921 got: $start, $end using $tempinfolog \n";

    $offset = $offset - $order;
    $end    = $end - $order;

    #    if ( $start < 1 ) {
    #        $start = 1;
    #    }
    $offset = $offset - $start + 1;
    print STDERR
"end:$end offset:$offset start:$start  donor:$donor accept:$accept ATG: $ATG\n";
    print STDERR "new offset: $offset\nnew start: $start\nnew order: $order\n";

    my $my_T_generic_lograt_summatrix_fn =
"${my_T_generic_logratio_freq_matrix_fn}_${order}_${matrix_type}.submatrix";
    my $my_T_generic_matrix_4param_fn =
"${my_T_generic_logratio_freq_matrix_fn}_${order}_${matrix_type}.matrix_4param";

#my $my_T_generic_lograt_summatrix_fn     = $my_T_logratio_freq_matrix_fn . ".submatrix";

#~ my $my_T_donor_dimatrix_4param_fn      = $my_T_logratio_freq_matrix_fn . ".dimatrix_4param";

#~ my $my_T_acceptor_lograt_summatrix_fn  = $my_T_logratio_freq_matrix_fn . ".submatrix";
#~ my $my_T_acceptor_dimatrix_4param_fn   = $my_T_logratio_freq_matrix_fn . ".dimatrix_4param";

#~ my $my_T_ATG_lograt_summatrix_fn       = $my_T_logratio_freq_matrix_fn . ".submatrix";
#~ my $my_T_ATG_dimatrix_4param_fn        = $my_T_logratio_freq_matrix_fn . ".dimatrix_4param";

    #~ my %my_dimatrics_dict = (
    #~ donor    => {pre => 2,
    #~ new => 3,
    #~ post => 4,
    #~ A_fn => $my_T_donor_lograt_summatrix_fn,
    #~ B_fn =>  $my_T_donor_dimatrix_4param_fn,
    #~ B_exec => $exec_B1 },
    #~ acceptor => 0.04,
    #~ ATG      => 0.15,
    #~ branch   => 0.30,
    #~ );

## DONOR DIMATRIX START
    if ( $order >= 1 && $donor ) {

        #my $pre_offset = $offset + $my_dimatrics_dict{'donor'};

        my $pre_offset  = $offset + 2;
        my $new_offset  = $offset + 3;
        my $post_offset = $offset + 4;

        #my_True_dimatrixdonor_4param_fn
        my $exec_A1 = "gawk -f ./bin/submatrix.awk ";
        my $exec_B1 = "./bin/preparedimatrixdonor4parameter.awk ";
        my $my_command_A =
"$exec_A1 $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_lograt_summatrix_fn";
        my $my_command_B =
"$exec_B1 $pre_offset $new_offset $post_offset $my_T_generic_lograt_summatrix_fn > $my_T_generic_matrix_4param_fn";

        print STDERR "$my_command_A \n";
        run($my_command_A);
        print STDERR "$my_command_B \n";
        run($my_command_B);

# print STDERR "submatrix.awk $start $end $true_seq_name-log.$ordname-matrix/$true_seq_name-log-info.$ordname-matrix";

    }
## DONOR DIMATRIX END
## ACCEPTOR DIMATRIX START
    elsif ( $order >= 1 && $accept ) {

        my $pre_offset  = $offset - 1;
        my $new_offset  = $offset;
        my $post_offset = $offset + 1;

        my $my_command_A =
"gawk -f ./bin/submatrix.awk $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_lograt_summatrix_fn";
        print STDERR "$my_command_A \n";
        run($my_command_A);

        my $my_command_B =
"./bin/preparedimatrixacceptor4parameter.awk $pre_offset $new_offset $post_offset $my_T_generic_lograt_summatrix_fn > $my_T_generic_matrix_4param_fn";
        print STDERR "$my_command_B \n";
        run($my_command_B);

#~ run(
#~ "./bin/preparedimatrixacceptor4parameter.awk $pre_offset $new_offset $post_offset $my_Tlograt_summatrix_fn > $my_T_dimatrixdonor_4param_fn"
#~ );

#~ #      print STDERR "submatrix.awk $start $end $true_seq_name-log.$ordname-matrix/$my_True_dimatrixdonor_4param_fn";

    }
## ACCEPTOR DIMATRIX END

## ATG DIMATRIX START
    elsif ( $order >= 2 && $ATG ) {

        my $pre_offset  = $offset - 2;
        my $new_offset  = $offset - 1;
        my $post_offset = $offset;

        run(
" gawk -f ./bin/submatrix.awk $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_lograt_summatrix_fn"
        );
        run(
" ./bin/preparetrimatrixstart4parameter.awk $pre_offset $new_offset $post_offset $my_T_generic_lograt_summatrix_fn > $my_T_generic_matrix_4param_fn"
        );
    }

## ATG DIMATRIX END

## ALL REMAINING CASES START
    else {

# print STDERR "$path/submatrix_order0.awk $start $end $true_seq_name-log.$ordname-matrix\n";

        run(
" gawk -f ./bin/submatrix_order0.awk $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_matrix_4param_fn"
        );

    }
## ALL REMAINING CASES END

## CREATE DATA STRUCTURE CONTAINING MATRIX OF INTEREST

    open( my $fh_PROF, "<", "$my_T_generic_matrix_4param_fn" )
      or croak "Failed here";
    while (<$fh_PROF>) {
        next if m/^#/;
        last if m/^\s/;
        last if m/^[^\d]/;
        chomp;
        my @e = split;
        push @prof, \@e;
    }
    close $fh_PROF;

    $prof_len = $end - $start + 1;
    print STDERR "length: $end - $start / $prof_len \n";
    return ( \@prof, $prof_len, $offset, $start, $end );

    #unlink $tempinfolog;
    #unlink "$my_T_dimatrixdonor_4param_fn";
    #unlink "$my_Tlograt_summatrix_fn";

}    ####END GETKMATRIX FUNCTION (Splice site an Start codon PWMs)

sub numerically { $a <=> $b }

sub fold4fasta {
    my $seq       = shift;
    my $foldedseq = "";

    #my $s = "";
    my $position_in_seq = 0;
    my $seq_len         = length($seq);

    while ( $position_in_seq < $seq_len ) {
        my $out_seq = substr( $seq, $position_in_seq, 60 );
        $foldedseq = $foldedseq . $out_seq . "\n";
        $position_in_seq += 60;
    }
    return $foldedseq;
}

#Optimize parameter file
sub OptimizeParameter {

    my (
        $gpfa,         $gpgff,        $newparam,      $branchswitch,
        $prof_len_bra, $fxdbraoffset, $branch_matrix, $IeWF,
        $deWF,         $FeWF,         $IoWF,          $doWF,
        $FoWF,         $iMin,         $dMin,          $fMin,
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

    print STDERR "\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\n\n";
    print {$fh_SOUT}
      "\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\n\n";

    for ( $IeWF = $IeWFini ; $IeWF <= $FeWF ; $IeWF += $deWF ) {    #for_#1
        print STDERR "eWF: $IeWF\noWF: ";

        for ( $IoWF = $IoWFini ; $IoWF <= $FoWF ; $IoWF += $doWF ) {    #for_#2
            print STDERR "$IoWF  ";
            my $param = Geneid::Param->new();
            $param->readParam("$newparam");

            for ( my $i = 0 ; $i < $param->numIsocores ; $i++ ) {       #for_#3
                if ( !defined @{ $param->isocores }[$i]
                    ->Exon_weights( [ $IeWF, $IeWF, $IeWF, $IeWF ] ) )
                {
                    croak "error in setting exon weights L3163\n";
                }

                if ( !defined @{ $param->isocores }[$i]
                    ->Exon_factor( [ $IoWF, $IoWF, $IoWF, $IoWF ] ) )
                {
                    croak "error in setting exon weights\n";
                }

#~ #   if (!defined @{$param->isocores}[$i]->Exon_factor([0.4,$IoWF,$IoWF,0.4])) {

                if (
                    !defined @{ $param->isocores }[$i]->Site_factor(
                        [ 1 - $IoWF, 1 - $IoWF, 1 - $IoWF, 1 - $IoWF ]
                    )
                  )
                {
                    croak "error in setting exon weights L3175 \n";
                }

#~ #      if (!defined @{$param->isocores}[$i]->Site_factor([0.55,1-$IoWF,1-$IoWF,0.55])) {
#~ }

            }    #end for_#3

            my $temp_geneid_param = "$work_dir/$species.geneid.param.tmp";
            $param->writeParam($temp_geneid_param);
            ###
            #~ my $fh_geneid    = File::Temp->new();
            #~ my $fname_geneid = $fh_geneid->filename;
            my $fh_geneid;
            my $fname_geneid;
            ( $fh_geneid, $fname_geneid ) =
              tempfile( DIR => $geneid_dir, SUFFIX => '.geneid' );
            print "\ntemp geneid file: $fname_geneid \n";
            my $my_command =
              "./bin/geneid -GP $temp_geneid_param $gpfa > $fname_geneid";
            run($my_command);
            my $temp0_geneid_pred_gff_fh =
              "$tmp_dir/Predictions." . basename($newparam) . ".gff";
            $my_command =
"cat $fname_geneid | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' > $temp0_geneid_pred_gff_fh ";
            run($my_command);

#` ./bin/geneid -GP ${newparam}.temp $gpfa | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' > $tmp_dir/Predictions.${newparam}.gff`;
## BUG very complex comand line

            my $temp0_evalout_fn =
              "$tmp_dir/Predictions." . basename($newparam) . ".tmp0_eval_out";
            $my_command =
"./bin/evaluation -sta $temp0_geneid_pred_gff_fh $gpgff  > $temp0_evalout_fn";
            print "\n$my_command\n";
            run($my_command);

            my $temp1_evalout_fn =
              "$tmp_dir/Predictions." . basename($newparam) . ".tmp1_eval_out";
            $my_command = "tail -2 $temp0_evalout_fn | head -1 |  
            gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $IoWF, $IeWF, \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}' > $temp1_evalout_fn ";
            print "\n$my_command\n";
            run($my_command);

            my @evaluation_output;
            open( my $fh_IN, "<", "$temp1_evalout_fn" ) or croak "Failed here";
            while (<$fh_IN>) {
                @evaluation_output = split " ";

                #@evaluation_output = split " ", $_;
            }
            close $fh_IN;

            #####

            #~ my @evaluation_output = split " ",

#~ #` ./bin/evaluation -sta $temp0_geneid_pred_gff_fh $gpgff | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $IoWF, $IeWF, \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}' `;

            push( @evaluation_total, \@evaluation_output );

        }    #end for_#2
        $IoWF = $IoWFini;
        print STDERR "\n";

    }    #end for_#1

    return \@evaluation_total;

}
## end sub optimize parameter file

sub BuildOptimizedParameterFile {

    my (
        $evalarray,    $branchswitch, $prof_len_bra,
        $fxdbraoffset, $branch_matrix
    ) = @_;
    my @sortedeval     = ();
    my @evaluationinit = ();
    my $best_IoWF      = "";
    my $best_IeWF      = "";
    my $best_Min       = "";
    my $best_Acc       = "";
    my $fh_SOUT;
    open( $fh_SOUT, ">", "$work_dir/$species.BuildOptimizedParameterFile.log" )
      or croak "Failed here";

    ## DEBUG
    print STDERR "\n input to OptimizeParameter sub
             \n1:"
      . $evalarray . "\n2 "
      . $branchswitch . "\n3 "
      . $prof_len_bra . "\n4 "
      . $fxdbraoffset . "\n5 "
      . $branch_matrix . "\n";

    if ( !$branchswitch ) {
        ## BUG ???
        @sortedeval = sort sorteval @{$evalarray};

        $best_IoWF = $sortedeval[0][0];    #0.2
        $best_IeWF = $sortedeval[0][1];    #-3.5

        print STDERR "\nBest performance obtained using IoWF: "
          . $sortedeval[0][0]
          . " and IeWF: "
          . $sortedeval[0][1] . "\n";
        print {$fh_SOUT} "\nBest parameter file performance obtained using oWF: "
          . $sortedeval[0][0]
          . " and eWF: "
          . $sortedeval[0][1] . "\n";

        #INITIALIZE ARRAY WITH EVALUATION PARAMETERS
        @evaluationinit =
          (qw(oWF eWF SN SP CC SNe SPe SNSP SNg SPg SNSPg raME raWE));

        print STDERR
"Sorted performance results (Three best performance estimates) for different values of oWF and eWF:\n"
          . join( "\t", @evaluationinit ), "\n";
        print $fh_SOUT
"Sorted performance results (best to worst) for different values of oWF and eWF:\n\n"
          . join( "\t", @evaluationinit ), "\n";

        foreach my $eval_ref (@sortedeval) {

            print $fh_SOUT join( "\t", @{$eval_ref} ), "\n";

        }

## FOUR BEST PERFORMANCE RESULTS AFTER OPTIMIZATION
        for ( my $i = 0 ; $i <= 2 ; $i++ ) {
            print STDERR join( "\t", @{ $sortedeval[$i] } ), "\n";
        }
############

## BUILD BEST PERFORMANCE (OPTIMIZED) PARAMETER FILE (USE VALUES OBTAINED ABOVE)

        my $param = Geneid::Param->new();
        $param->readParam("$newparam");

        for ( my $i = 0 ; $i < $param->numIsocores ; $i++ ) {
            if (
                !defined @{ $param->isocores }[$i]->Exon_weights(
                    [ $best_IeWF, $best_IeWF, $best_IeWF, $best_IeWF ]
                )
              )
            {
                croak "error in setting exon weights\n";
            }
            if (
                !defined @{ $param->isocores }[$i]->Exon_factor(
                    [ $best_IoWF, $best_IoWF, $best_IoWF, $best_IoWF ]
                )
              )
            {
#     if (!defined @{$param->isocores}[$i]->Exon_factor([0.4,$best_IoWF,$best_IoWF,0.4])) {
                croak "error in setting exon weights\n";
            }
            if (
                !defined @{ $param->isocores }[$i]->Site_factor(
                    [
                        1 - $best_IoWF,
                        1 - $best_IoWF,
                        1 - $best_IoWF,
                        1 - $best_IoWF
                    ]
                )
              )
            {
#     if (!defined @{$param->isocores}[$i]->Site_factor([0.55,1-$best_IoWF,1-$best_IoWF,0.55])) {
                croak "error in setting exon weights\n";
            }
        }

        #write new parameter file (optimized)
        my $optimized_geneid_param_fn =
          "$results_dir/$species.geneid.optimized.param";
        $param->writeParam("$species.geneid.optimized.param");

        print STDERR
"\n\nNew optimized parameter file named: $species.geneid.optimized.param \n";
        print $fh_SOUT
"\n\nNew optimized parameter file named: $species.geneid.optimized.param \n";

        return [ $best_IeWF, $best_IoWF, 0, 0, \@evaluationinit ];

    }
    close $fh_SOUT;
    return 1;
}

sub EvaluateParameter {

    my ( $gpfa, $gpgff, $newparam, $IoWF, $IeWF ) = @_;
    my $my_command;

    my $geneid_test_predict_gff_fn =
      "$work_dir/geneid_test_predictions." . basename($newparam) . ".gff";
    print STDERR "\ngeneid_test_predict_gff_fn: $geneid_test_predict_gff_fn\n";

    #~ my $fh_geneid;
    #~ my $fname_geneid;
    my ( $fh_geneid, $fname_geneid ) =
      tempfile( DIR => $geneid_dir, SUFFIX => '.eval_par.geneid' );
    print "\ntemp geneid file: $fname_geneid \n";
    $my_command = "./bin/geneid -GP $newparam $gpfa > $fname_geneid";

    #print STDERR "\n$my_command, not running\n";
    run($my_command);

    #` ./bin/geneid -GP $newparam $gpfa
    $my_command =
"cat $fname_geneid | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' > $geneid_test_predict_gff_fn";
    run($my_command);

###dk
#` ./bin/geneid -GP ${newparam}.temp $gpfa | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' > $tmp_dir/Predictions.${newparam}.gff`;
## BUG very complex comand line

    my $tempA_evalout_fn =
      "$geneid_dir/Predictions." . basename($newparam) . ".tmpA_eval_out";
    $my_command =
"./bin/evaluation -sta $geneid_test_predict_gff_fn $gpgff > $tempA_evalout_fn";
    print "\n$my_command\n";
    run($my_command);

    my $tempB_evalout_fn =
      "$geneid_dir/Predictions." . basename($newparam) . ".tmpB_eval_out";
    $my_command = "tail -2 $tempA_evalout_fn | head -1 |  
            gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $IoWF, $IeWF, \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}' > $tempB_evalout_fn ";
    print "\n$my_command\n";
    run($my_command);

    my @evaluation_test;
    open( my $fh_IN, "<", "$tempB_evalout_fn" ) or croak "Failed here";
    while (<$fh_IN>) {
        @evaluation_test = split " ";

        #@evaluation_output = split " ", $_;
    }
    close $fh_IN;
###dk

#~ my @evaluation_test = split " ",
#~ ` ./bin/evaluation -sta $geneid_test_predict_gff_fn $gpgff | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}' `;

    return \@evaluation_test;

}    # evaluate parameter function

sub calculate_stats {

    my (
        $species,                  $sout,
        $out_intron_X,             $cds_all_nozero_tbl,
        $out_gff_X,                $inframe_X,
        $inframe_X_eval,           $seqs_used_XX,
        $tot_noncanon_donors_intX, $tot_noncanon_accept_intX,
        $tot_noncanon_ATGx,        $markov_model,
        $total_coding,             $total_noncoding,
        $stdo,                     $endo,
        $stac,                     $enac,
        $stst,                     $enst,
        $stbr,                     $enbr,
        $branchsw,                 $use_allseqs_flag
    ) = @_;
    my $my_command = "";

## OBTAIN GENE MODEL SET STATISTICS
## Open gene model object

    $param->geneModel( Geneid::GeneModel->new() );
    $param->geneModel->useDefault;
###########
    my $fh_SOUT;
    open( $fh_SOUT, ">", "$work_dir/test.WriteStatsFileLog.txt" )
      or croak "Failed here";

    my $avgintron = "";
    my $sdintron  = "";

    #my @intronlength = ` gawk '{print length(\$2)}' $out_intron_X `;

    #print STDERR "INTRON: $mean, $st\n";
    $my_command = "gawk '{print length(\$2)}' $out_intron_X | sort -n";
    my @introns_len_list = capture($my_command);
    my ( $mean, $st ) = average( \@introns_len_list );
    print STDERR "INTRONS mean, ST: $mean, $st\n";

    #~ ## BUG wrong command?

    #~ #$my_command = "gawk '{print length(\$2)}' $out_intron_X | sort -n";
    #~ $my_command = "sort -k2,2n $out_intron_X";
    #~ my @intronlist = capture($my_command);

    my $totintrons = scalar(@introns_len_list);
    my @intronlen  = ();
    my $intr       = "";
    for ( my $i = 0 ; $i <= scalar(@introns_len_list) - 1 ; $i++ ) {

        $intr = $introns_len_list[$i];
        chomp $intr;
        push( @intronlen, $intr );
    }

    my @slice1 = @intronlen[ 0 .. 5 ];
    my @slice2 =
      @intronlen[ ( scalar(@intronlen) - 5 ) .. ( scalar(@intronlen) - 1 ) ];

    ## BUG 2857
    my $intron_short_int = $introns_len_list[0] - $introns_len_list[0] * (0.25);
    chomp $intron_short_int;
    if ( $intron_short_int > 40 ) { $intron_short_int = 40; }

    #my $intron_long_int =  $intronlist[$totintrons - 1];
    my $intron_long_int =
      $mean + ( $st * 3 ) > 100_000 ? 100_000 : $mean + ( $st * 3 );
    chomp $intron_long_int;

    my $intergenic_min = 200;
    my $intergenic_max = 'Infinity';

## use shortest and longest intron lengths in gene model of parameter file
    $param->geneModel->intronRange( $intron_short_int, $intron_long_int );
    $param->geneModel->intergenicRange( $intergenic_min, $intergenic_max );
###############################

    $my_command =
      "gawk '{print gsub(/[GC]/,\".\",\$2)/length(\$2)}' $cds_all_nozero_tbl";
    my @CDSGCcontent = capture($my_command);

    #print STDERR "@CDSGCcontent\n";
    my ( $meangc, $stgc ) = average( \@CDSGCcontent );

    #print STDERR "CDS: $meangc $stgc $out_intron_X\n";
    $my_command =
      "gawk '{print gsub(/[GC]/,\".\",\$2)/length(\$2)}' $out_intron_X ";
    my @intronGCcontent = capture($my_command);

    #print STDERR "@intronGCcontent\n";
    my ( $meangci, $stgci ) = average( \@intronGCcontent );

 #print STDERR "intron: $meangci $stgci\n";
 #BUG?
 #my $totexons = ` gawk '{print \$9}' $out_gff_X | wc -l | gawk '{print \$1}' `;
    $my_command = "gawk '{print \$9}' $out_gff_X | wc -l ";
    my $totexons = capture($my_command);

    #my $totexons = ` gawk '{print \$9}' $out_gff_X | wc -l `;

    chomp $totexons;
    $totexons = int($totexons);
    my @exonspergene;
    $my_command =
      "gawk '{print \$9}' $out_gff_X | sort | uniq -c | gawk '{print \$1}'";
    @exonspergene = capture($my_command);

    #@exonspergene =
    #  ` gawk '{print \$9}' $out_gff_X | sort | uniq -c | gawk '{print \$1}' `;

    my ( $avgex, $stex ) = average( \@exonspergene );

    $my_command =
      "egrep -v 'Single' $out_gff_X | gawk '{len=\$5-\$4;print len}' - | sort ";
    my @exonlength = capture($my_command);

#    my @exonlength =
#      ` egrep -v 'Single' $out_gff_X | gawk '{len=\$5-\$4;print len}' - | sort `;

    my ( $avgle, $stle ) = average( \@exonlength );

    $my_command = "egrep -c '(Single)' $out_gff_X";
    my $singlegenes = capture($my_command);

    #my $singlegenes = `egrep -c '(Single)' $out_gff_X `;
    chomp $singlegenes;

    #print $fh_SOUT "GENE MODEL STATISTICS FOR $species\n\n";

    print $fh_SOUT
"\nA subset of $tot_seqs4training_intX sequences (randomly chosen from the $total_seqs gene models) was used for training\n\n";
    my $total_seqs;
    if ( !$use_allseqs_flag ) {
        print $fh_SOUT
"The user has selected to use $seqs_used_XX gene models (80 % of total) for training and to set aside BUG  was xxgffseqseval annotations (20 % of total) for evaluation\n\n";
    }
    else {
        print $fh_SOUT
"$total_seqs gene models were used for both training and evaluation\n\n";
    }

    if ( !$use_allseqs_flag ) {
        print $fh_SOUT
"$inframe_X of the gene models translate into proteins with in-frame stops within the training set and $inframe_X_eval in the evaluation set (seqs removed).\n\n";
    }
    else {
        print $fh_SOUT
"$inframe_X of the gene models translate into proteins with in-frame stops within the training set.\n\n";
    }
    print $fh_SOUT
"There are $tot_noncanon_donors_intX non-canonical donors as part of the training set\n\n";
    print $fh_SOUT
"There are $tot_noncanon_accept_intX non-canonical acceptors as part of the training set\n\n";
    print $fh_SOUT
"There are $tot_noncanon_ATGx non-canonical start sites as part of the training set\n\n";
    print $fh_SOUT
"These gene models correspond to $total_coding coding bases and $total_noncoding non-coding bases\n\n";
    print $fh_SOUT
"Deriving a markov model for the coding potential of order $markov_model\n\n";
    print $fh_SOUT
"The intronic sequences extracted from the gene models have an average length of $mean, with $st of SD\n";
    print $fh_SOUT
"Geneid can predict gene models having introns with a minimum length of $intron_short_int nucleotides and a maximum of $intron_long_int bases (boundaries used in gene model) \n\n";
    print $fh_SOUT
"The minimum (user selected) intergenic distance was set to $intergenic_min nucleotides whereas the maximum was set to $intergenic_max (boundaries used in gene model) \n\n";
    print $fh_SOUT
"The GC content of the exonic and intronic sequences is $meangc (SD $stgc) and $meangci (SD $stgci) respectively \n\n";
    print $fh_SOUT
      "The gene models used for training contain $totexons exons \n\n";
    print $fh_SOUT
      "The gene models average $avgex exons per gene (SD $stex)\n\n";
    print $fh_SOUT
"The average length of the exons (non-single) in the training set gene models is $avgle (SD $stle)\n\n";

    if ( !$use_allseqs_flag ) {
        print $fh_SOUT
"The training set includes $singlegenes single-exon genes (out of $seqs_used_XX ) gene models\n\n";
    }
    else {
        print $fh_SOUT
"The training set includes $singlegenes single-exon genes (out of $total_seqs) gene models\n\n";
    }
    print $fh_SOUT "The donor site profile chosen by the user spans "
      . ( $endo - $stdo + 1 )
      . " nucleotides: position $stdo to $endo\n";
    print $fh_SOUT "The acceptor site profile chosen by the user spans "
      . ( $enac - $stac + 1 )
      . " nucleotides: position $stac to $enac\n";
    print $fh_SOUT "The start site profile chosen by the user spans "
      . ( $enst - $stst + 1 )
      . " nucleotides: position $stst to $enst\n";

    if ($branchsw) {
        print $fh_SOUT "The branch site profile chosen by the user spans "
          . ( $enbr - $stbr + 1 )
          . " nucleotides: position $stbr to $enbr\n";
    }
    close $fh_SOUT;
    return ( $intron_short_int, $intron_long_int, $intergenic_min,
        $intergenic_max );

}

sub average {
    my ($sequences) = @_;
    my $sum         = 0;
    my $total       = 0;
    my ( $mean, $st );

    foreach my $seq ( @{$sequences} ) {
        $sum += $seq;
        $total++;
    }

    $mean = $sum / $total;
    $mean = sprintf( "%.3f", $mean );

    $sum = 0;

    foreach my $seq ( @{$sequences} ) {
        $sum += ( $seq - $mean ) * ( $seq - $mean );

    }
    $st = sqrt( $sum / $total );
    $st = sprintf( "%.3f", $st );

    return ( $mean, $st );
}

sub TblToFasta {
    my ( $tbl, $faout ) = @_;

    open( my $fh_IN,   "<", "$tbl" )   or croak "Failed here";
    open( my $fh_FOUT, ">", "$faout" ) or croak "Failed here";
    while (<$fh_IN>) {
        chomp;
        #~ my ( $n, $s ) = split( /\s+/, $_ );
        my ( $n, $s ) = split( /\s+/);
        my ( $i, $e ) = ( 1, length($s) );
        print {$fh_FOUT} ">$n\n";
        while ( $i <= $e ) {
            print {$fh_FOUT} substr( $s, $i - 1, 60 ) . "\n";
            $i += 60;
        }
    }
    close $fh_IN;
    close $fh_FOUT;
    return $faout;
}

sub TblToFastaFile {

    my ( $dir, $tbl_fn ) = @_;

    open( my $fh_IN_tbl, "<", "$tbl_fn" ) or croak "Failed here";

    print STDERR "## $tbl_fn\n";
    while (<$fh_IN_tbl>) {
        my $input_line = $_;
        chomp($input_line);

        #my ( $seq_name, $seq ) = split( /\s+/o, $_ );
        my ( $seq_name, $seq ) = split( /\s+/o, $input_line );

        #print STDERR "YYY $seq_name \t";
        ##open( FOUT, ">${dir}" . "$n" );
        open( my $fh_FOUT_fasta, ">", "${dir}" . "$seq_name" );

        #print STDERR "XXX ${dir}" . "$seq_name \n";
        my ( $base_number, $seq_length ) = ( 1, length($seq) );
        print {$fh_FOUT_fasta} ">$seq_name\n";
        print STDERR "#";
        while ( $base_number <= $seq_length ) {
            print {$fh_FOUT_fasta} substr( $seq, $base_number - 1, 60 ) . "\n";
            $base_number += 60;
        }
        close $fh_FOUT_fasta;
    }

    close $fh_IN_tbl;

    #close $fh_FOUT_fasta;
    return 1;
}

sub fasta_2_tbl {
    my ( $fa, $tblout ) = @_;

    open( my $fh_IN,   "<", "$fa" );
    open( my $fh_TOUT, ">", "$tblout" );

    #print STDERR "$fa INSIDE LOOP\n\n";
    my $count = 0;
    while (<$fh_IN>) {
        chomp;
        $_ =~ s/\|//;
        if ( $_ =~ /\>(\S+)/ ) {
            if ( $count > 0 ) {
                print {$fh_TOUT} "\n";
            }
            print {$fh_TOUT} $1 . "\t";
            $count++;
        }
        else {
            print {$fh_TOUT} $_;
        }
    }
    print {$fh_TOUT} "\n";

    close $fh_IN;
    close $fh_TOUT;

    return $tblout;
}

#~ sub FastaToTbl {

#~ my ( $in_fa_fn, $out_tbl_fn, $flag ) = @_;
#~ say "\n$in_fa_fn, $out_tbl_fn, $flag\n";
#~ #fix CAPS & Ns while converting
#~ #

#~ open( my $fh_IN,   "<", "$in_fa_fn" ) ;
#~ open( my $fh_TOUT, ">", "$out_tbl_fn" );

#~ #print STDERR "$fa INSIDE LOOP\n\n";
#~ my $count    = 0;
#~ my $sequence = "";
#~ while (<$fh_IN>) {
#~ chomp;
#~ $_ =~ s/\|//;
#~ if ( $_ =~ /\>(\S+)/ ) {
#~ print $fh_TOUT "\n" if $count > 0;
#~ print $fh_TOUT $1 . "\t";
#~ $count++;
#~ }
#~ else {
#~ if ($flag eq "genome"){
#~ say "here\n";
#~ $sequence = uc $_;
#~ #$sequence = ~ s/-/N/gi;
#~ print $fh_TOUT $_;
#~ }
#~ else {
#~ print $fh_TOUT $_;
#~ }
#~ }
#~ }
#~ print $fh_TOUT "\n";

#~ close $fh_IN;
#~ close $fh_TOUT;
#~ if ($flag eq "foobar"){
#~ my $my_command = "sort --output=$out_tbl_fn $out_tbl_fn";
#~ run($my_command);
#~ }
#~ #return $tblout;
#~ return 1;

#~ }

sub translate_2_protein {

    my ( $geneticcode_fn, $cds_fn, $outprot ) = @_;
    print STDERR "$geneticcode_fn in loop\n";
    my $frame = 0;

    #~ if ( !open(my $fh_FILEIN, "<", "$geneticcode" ) ) {
    #~ print "translate_2_protein: impossible to open genetic.code\n";
    #~ exit(1);
    #~ }

    my %gencodeh = ();
    open( my $fh_gencode, "<", "$geneticcode_fn" )
      or croak "Can't open  $geneticcode_fn";
    while (<$fh_gencode>) {

        my $line = $_;

        my ( $aa, $codon ) = split( /\s+/, $line );

        #print STDERR "$codon\n";

        $gencodeh{$codon} = $aa;

        #print STDERR "HERE: $gencodeh{$codon}\n";
    }
    close $fh_gencode;

    #if ( !open(my $fh_CDSIN, "<", "$cds" ) ) {
    #print "translate_2_protein: impossible to open $cds\n";
    #exit(1);
    #}

    open( my $fh_CDS_IN, "<", "$cds_fn" )  or croak "err:  $cds_fn";
    open( my $fh_POUT,   ">", "$outprot" ) or croak "Failed here";
    print STDERR "translating: $cds_fn \n";

    while (<$fh_CDS_IN>) {

        my $line = $_;

        my ( $name, $seq ) = split( /\s+/, $line );

        my $lseq = length($seq);

        print {$fh_POUT} "$name ";

        for ( my $i = $frame ; $i < $lseq - 1 ; $i += 3 ) {
            my $triplet = substr( $seq, $i, 3 );
            if ( $gencodeh{$triplet} ) {
                print {$fh_POUT} "$gencodeh{$triplet}";
            }
            else {
                print {$fh_POUT} "X";
            }

        }    #for
        print {$fh_POUT} "\n";
    }

    close $fh_CDS_IN;
    close $fh_POUT;
    say "\noutprot : $outprot \n";
    return $outprot;

}

## CONVERT GFF2 TO GENEID GFF
sub convert_GFF_2_geneidGFF {

    my ( $no_dots_gff_fn, $species, $type ) = @_;
    my %G;
    my @G = ();

    my $geneid_gff = $work_dir . $species . ${type} . ".geneid_gff";

    open( my $fh_GFF,    "<", "$no_dots_gff_fn" ) or croak "Failed here";
    open( my $fh_GFFOUT, ">", "$geneid_gff" )      or croak "Failed here";
    while (<$fh_GFF>) {
        my ( $c, @f, $id );
        $c = ":";
        $_ =~ s/\|//;
        chomp;
        @f = split(/\s+/o);
        $id = $f[8];    #seq name i.e. 7000000188934730
        ( exists( $G{$id} ) ) || do {
            $c = "#";
            $G{$id} = [ @f[ 8, 0, 6 ], 0, @f[ 3, 4 ], [] ]
              ;    # [7000000188934730 1.4 - 0 46549    46680 ]
        };
        push @{ $G{$id}[6] }, [ @f[ 3, 4 ] ];
        $G{$id}[3]++;
        $G{$id}[4] > $f[3] && ( $G{$id}[4] = $f[3] );
        $G{$id}[5] < $f[4] && ( $G{$id}[5] = $f[4] );
    }

    @G =
      sort { $a->[1] cmp $b->[1] || $a->[4] <=> $b->[4] || $a->[5] <=> $b->[5] }
      map { $G{$_} } keys %G;

    foreach my $GN (@G) {
        my (
            $id,     $ctg, $str, $nex, $go,   $ge, $coords,
            @coords, $ce,  $CE,  $cur, $feat, $c
        );
        ( $id, $ctg, $str, $nex, $go, $ge, $coords ) = @{$GN};

       # print STDERR Data::Dumper->Dump([ \@$coords ], [ qw/ *coords / ])."\n";
        @coords = map { $_->[0], $_->[1] }
          sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
          map { $_ } @{$coords};

        # print STDERR Data::Dumper->Dump([ \@coords ], [ qw/ *Coords / ])."\n";
        #print "# $id $ctg $str $nex $go $ge\n";
        $ce = 0;
        $CE = $nex * 2;
        $c  = 1;
        while ( $ce < $CE ) {

            # $cur = ($str eq '-' ? $CE - $ce - 2 : $ce);
            if ( $nex == 1 ) {
                $feat = "Single";
            }
            elsif ( $c == 1 ) {
                $feat = $str eq '-' ? "Terminal" : "First";
            }
            elsif ( $c == $nex ) {
                $feat = $str eq '-' ? "First" : "Terminal";
            }
            else {
                $feat = "Internal";
            }
            print {$fh_GFFOUT} join( "\t",
                $ctg, "$species", $feat, $coords[$ce], $coords[ ( $ce + 1 ) ],
                ".", $str, ".", $id )
              . "\n";
            $ce += 2;
            $c++;
        }
    }
    close $fh_GFF;
    close $fh_GFFOUT;

#~ my $X_geneid_sorted_gff = "";
#~ open( my $fh_LOCID, "-|", "sort -s -k8,9 -k4,5n $geneid_gff") or croak "Failed here";
#~ while (<$fh_LOCID>) {
#~ $X_geneid_sorted_gff .= $_;
#~ }
#~ close $fh_LOCID;

####
    my $geneid_sorted_gff_fn =
      $work_dir . $species . $type . ".geneid.gff_sorted";
    my $my_command = "sort -s -k8,9 -k4,5n $geneid_gff > $geneid_sorted_gff_fn";
    run($my_command);

    #~ open( my $fh_FOUT, ">", "$geneid_sorted_gff_fn" ) or croak "Failed here";
    #~ print $fh_FOUT "$geneid_sorted_gff_fn";
    #~ close $fh_FOUT;

    return $geneid_sorted_gff_fn;

    #exit(0);

}

sub sorteval {

         $b->[7] <=> $a->[7]
      || $b->[10] <=> $a->[10]
      || $b->[4] <=> $a->[4]
      || $a->[11] <=> $b->[11]
      || $a->[12] <=> $b->[12]

}

#DK_subs
sub num_of_lines_in_file {
    my $input_fn = $_[0];
    my $my_num_lines;
    $my_num_lines = capture("cat $input_fn | wc -l");

    #assert No such file or directory #
    chomp $my_num_lines;
    $my_num_lines = int($my_num_lines);
    return $my_num_lines;
}

sub check_external_progs() {
## Checking if the external programs are in the path.
## C, awk, python programs
    #my $prog_name;
    my $my_command;
    my @progs_2_check = (
        qw(bash
          gawk
          egrep
          sed
          geneid
          ssgff
          shuf
          pictogram
          gff2gp.awk
          cds2gff.awk
          frequency.py
          information.py
          submatrix.awk
          submatrix_order0.awk
          Getkmatrix.awk
          multiple_annot2one.awk
          logratio_kmatrix.awk
          logratio_zero_order.awk
          preparedimatrixacceptor4parameter.awk
          preparedimatrixdonor4parameter.awk
          preparetrimatrixstart4parameter.awk)
    );

    for my $prog_name (@progs_2_check) {
        $my_command = "which $prog_name > /dev/null";
        run($my_command);

        #or croak "\n $prog_name not available\n";
    }

    #BASH AWK
    #system("which gff2ps > /dev/null;")
    #    && &go_to_die("The gff2ps package is not found or is not executable");
    #say "\nNeccessary binaries are executable\n";
    return 1;
}

sub create_data_dirs {

    #my $dir_name = shift(@_);
    say "\ninside create_data_dirs function\n";
    ## check if not a bug XXX
    foreach my $dir_name (@_) {
        say "Creating $dir_name directory\n";
        if ( -d $dir_name ) {
            say "$dir_name exists. Purge and recreate\n";
            rmtree( ["dir_name"] );
        }
        else {
            #say "$dir_name does not exist!\n";
            my $my_command = "mkdir -p $dir_name";
            run($my_command);
        }
    }
    return 1;
}

sub write_sizes_from_tbl_fn {

    my $input_tbl_fn = $_[0];
    print "calc size for  $input_tbl_fn \n";

    open( my $fh_input, '<', $input_tbl_fn )
      or croak "Could not open file $input_tbl_fn' $OS_ERROR";

    while ( my $line = <$fh_input> ) {
        chomp $line;

        #print $line
        my @f = split / /, $line;
        my $name = $f[0];

        #print "$name\t"
        my $seq       = $f[1];
        my $lengfasta = length($seq);

        #print "$name = $lengfasta\t"
        my $tmp_out_fn = "$fastas_dir/$name" . "_len";

        #print  " $tmp_out_fn \t"
        open( my $fh_out, '>', $tmp_out_fn )
          or croak "Could not open file '$tmp_out_fn' $OS_ERROR";
        print {$fh_out} "$name $lengfasta\n";
        close $fh_out;

    }
    close $fh_input;
    return 1;
}    #end write_sizes_from_tbl_fn

sub get_background_kmers {

    my ( $kmer_len, $input_fas_fn, $tbl, $num_seqs, $background_tbl ) = @_;
    my $bckgrnd_debug_flag;
    $bckgrnd_debug_flag = 1;
    if ( !$bckgrnd_debug_flag ) {

        #my $input_fas_fn = $f;
        #my $tbl = "";
        my $countlines = 0;
        my $total_seqs = num_of_lines_in_file($input_fas_fn);

        #my $totalseqs  = `egrep -c \"^>\" $input_fas_fn`;
        #chomp $totalseqs;

        open( my $fh_IN, "<", "$tbl" ) or croak "Failed here";
        my @tabular = ();
        while (<$fh_IN>) {
            push @tabular, $_;
        }
        close $fh_IN;

        print STDERR
"\nThe total number of genomic sequences (in $input_fas_fn) is $total_seqs\n";

        if ( $total_seqs >= 1 ) {

            #print STDERR "in totalseqs if";
            my $seq    = "";
            my $sublen = 0;
            foreach my $line (@tabular) {
                chomp $line;
                my @f = split " ", $line;
                $seq .= $f[1];
                $sublen += length( $f[1] );
                $countlines++;
                if ( $sublen >= ( $num_seqs * $kmer_len + 1 ) ) {
                    last;
                }
            }

            #print STDERR "$seq";
            my $len = length($seq);

            # chomp $totalseqs;
            print STDERR
"The length of the new fasta file is $len\n(concatenation of $countlines fastas (out of $total_seqs))\n";
            open( my $fh_BACKGRND, ">", "$background_tbl" )
              or croak "Failed here";
            my $row = 1;
            print STDERR
"\nObtain background ($num_seqs seqs) from fasta of $kmer_len nucleotides \n";
            for ( my $n = 1 ; $n <= $num_seqs ; $n++ ) {
                my $kmer = "N";
                while ( $kmer =~ m/[^ACGTacgt]/ ) {
                    my $rand = int( rand( $len - $kmer_len ) + 1 );
                    $kmer = substr( $seq, $rand, $kmer_len );
                }

                #print STDERR  $n."..";
                print {$fh_BACKGRND} $row++ . "\t$kmer\n";
            }
            print STDERR "FINISHED OBTAINING BACKGROUND INFO\n";
            close $fh_BACKGRND;
        }
    }
    else {
        #just for non-random results
        print "\npremade background\n";
        my $my_command =
          "zcat ./test_data/pmar_1M_60mers.tbl.gz > $background_tbl";
        run($my_command);
    }

    return 1;

}    # END get_background_kmers function

sub compute_sites_pictogram {
#############################################################

    my $my_command;
    my $fh_FOUT;
#########
## get donor site statistics
#########
    my $order = "0";
    my $numbersites;
    $numbersites = num_of_lines_in_file($out_donor_tbl);

    my $don_offset =
      $bases_offset;   #position before intron (last of exon (31) -1 for offset)

    if ( $numbersites > $train_sites_cutoff ) {
        $order = "1";
    }
    elsif ( $numbersites <= $train_sites_cutoff ) {
        $order = "0";
    }

    print STDERR
"\nThere are $numbersites donor sites, enough for a matrix of order $order, prior offset: $don_offset $out_donor_tbl \n";

    my (
        $donor_matrix, $prof_len_don, $fxd_don_offset,
        $donor_start,  $donor_end
      )
      = get_K_matrix( $out_donor_tbl, $backgrnd_kmers_fn, $order, $don_offset,
        1, 0, 0, 0, 0, 0, 0 );
    if (
        !defined @{ $param->isocores }[0]->set_profile(
            'Donor_profile', $prof_len_don, $fxd_don_offset, $pwm_cutoff,
            $order, 0, 1, 0, 0, 0, 0, $donor_matrix
        )
      )
    {
        croak "error in setting profile\n";
    }

    my $donsub = "";

#print STDERR "gawk '{print  substr(\$2,($donor_start-3),($prof_len_don+6))}' $out_donor_tbl\n";

    $my_command =
"gawk '{print  substr(\$2,($donor_start-3),($prof_len_don+6))}' $out_donor_tbl ";
    $donsub = capture($my_command);

    $subprofile_donors = $work_dir . $species . ".don.sub.profile";
    open( $fh_FOUT, ">", "$subprofile_donors" ) or croak "Failed here";
    print {$fh_FOUT} "$donsub";
    close $fh_FOUT;

    $my_command =
"./bin/pictogram $subprofile_donors $plots_dir/donor_profile.pictogram -bits -land";
    print "\n$my_command\n";
    run($my_command);

#########
## get acceptor site statistics
#########

    $order = "0";

    $numbersites = num_of_lines_in_file($out_acceptor_tbl);

    #$numbersites = `wc -l $out_acceptor_tbl | gawk '{print \$1}'`;
    #chomp $numbersites;
    #$numbersites = int($numbersites);
    print "numbersites in $out_acceptor_tbl: $numbersites\n";
    my $acc_offset =
      $bases_offset;   #position after intron (first of exon (31) -1 for offset)

    if ( $numbersites > $train_sites_cutoff ) {

        $order = "1";

    }
    elsif ( $numbersites <= $train_sites_cutoff ) {

        $order = "0";
    }

    print STDERR
"\nThere are $numbersites acceptor sites, enough for a matrix of order $order, offset: $acc_offset \n";

    my (
        $acceptor_matrix, $prof_len_acc, $fxd_acc_offset,
        $acceptor_start,  $acceptor_end
      )
      = get_K_matrix( $out_acceptor_tbl, $backgrnd_kmers_fn, $order,
        $acc_offset, 0, 1, 0, 0, 0, 0, 0 );
    if (
        !defined @{ $param->isocores }[0]->set_profile(
            'Acceptor_profile', $prof_len_acc, $fxd_acc_offset, $pwm_cutoff,
            $order, 0, 1, 0, 0, 0, 0, $acceptor_matrix
        )
      )
    {
        croak "error in setting profile\n";
    }

    my $accsub = "";

    $my_command =
"gawk '{print  substr(\$2,($acceptor_start-3),($prof_len_acc+6))}' $out_acceptor_tbl ";
    $accsub = capture($my_command);

    $subprofile_acceptors = $work_dir . $species . ".acc.sub.profile";
    open( $fh_FOUT, ">", "$subprofile_acceptors" ) or croak "Failed here";
    print {$fh_FOUT} "$accsub";
    close $fh_FOUT;

    $my_command =
"./bin/pictogram $subprofile_acceptors $plots_dir/acceptor_profile.pictogram -bits -land";
    print "\n$my_command\n";
    run($my_command);

#########
## get start site statistics
#########

    $order       = "0";
    $numbersites = num_of_lines_in_file($out_ATGxs_tbl);

    my $ATGx_offset =
      $bases_offset;  #before first position of the exon (31)minus 1 for offset)

    if ( $numbersites > $train_sites_markov_cutoff ) {

        $order = "2";

    }
    elsif ( $numbersites <= $train_sites_markov_cutoff ) {

        $order = "0";
    }

    print STDERR
"\nThere are $numbersites ATGx start sites, enough for a matrix of order $order, offset: $ATGx_offset \n";

    my ( $start_matrix, $prof_len_sta, $fxd_ATGx_offset, $ATGx_start,
        $ATGx_end ) =
      get_K_matrix( $out_ATGxs_tbl, $backgrnd_kmers_fn, $order, $ATGx_offset,
        0, 0, 1, 0, 0, 0, 0 );

## write to parameter file
    if (
        !defined @{ $param->isocores }[0]->set_profile(
            'Start_profile', $prof_len_sta, $fxd_ATGx_offset, $pwm_cutoff,
            $order, 0, 1, 0, 0, 0, 0, $start_matrix
        )
      )
    {
        croak "error in setting profile\n";
    }
#############################

    my $stasub = "";

    $my_command =
"gawk '{print  substr(\$2,($ATGx_start-3),($prof_len_sta+6))}' $out_ATGxs_tbl ";
    $stasub = capture($my_command);

    $subprofile_ATGs = $work_dir . $species . ".ATGx.sub.profile";
    open( $fh_FOUT, ">", "$subprofile_ATGs" ) or croak "Failed here";
    print {$fh_FOUT} "$stasub";
    close $fh_FOUT;

    $my_command =
"./bin/pictogram $subprofile_ATGs $plots_dir/ATGx_profile.pictogram -bits -land";

    print "\n$my_command\n";
    run($my_command);

## end get start site statistics
    return $donor_start, $donor_end, $acceptor_start, $acceptor_end,
      $ATGx_start,
      $ATGx_end;

#############################################################
}
