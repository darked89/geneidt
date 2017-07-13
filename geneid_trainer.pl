#!/usr/bin/perl

## if run under perlbrew, use i.e.:
#!/$HOME/perl5/perlbrew/perls/perl-5.10.1/bin/perl

## checks & debugs modules
use Modern::Perl;
use English '-no_match_vars';

## use experimental 'signatures';

no warnings "experimental::signatures";
use feature "signatures";

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
## 2016.12.11 use Benchmark qw(:all);
##use YAML qw(Dump Bless);
use YAML::XS qw(LoadFile);
use Smart::Comments;

##use Function::Parameters qw(:strict);
## Problem with my 5.18.1 install @CRG
## use Devel::Size qw(size total_size)

## geneid_trained modules
use lib '.';
use Geneid::Param;
use Geneid::Isocore;
use Geneid::geneid;
use Geneid::geneidCEGMA;

## experimental
## use Inline::Python;

## MAIN VARIABLES
## my $PROGRAM_NAME    = "geneid_trainer";
## my $PROGRAM_VERSION = "2016.12.15a";
my $PROGRAM_HOME = getcwd;

my $exec_path = "$PROGRAM_HOME/bin/";

local $ENV;
$ENV{'PATH'} = $exec_path . ":" . $ENV{'PATH'};

## TODO: use YAML to get results from Python pre-processor script
#~ my $config = LoadFile('./dev/conf_001.yaml');
#~ print Dumper($config);
#~ say $config->{genome_fasta};
#~ die("checking YAML");

my $genetic_code = "./etc/genetic.code";

### $genetic_code geneticcode at <file>[<line>]

## no need to run anything if this fails
check_external_progs();

## Move parts necessary for getting comand line args here
my $species      = "";
my $input_gff_fn = "";
my $input_fas_fn = "";
my $sout         = "-";
my $fix_fasta    = 0;

## my $branch_pred_flag = 0;

## my $reduced     = 0;
## my $interactive_flag = 0;

#~ my $tenfold     = 0;
#~ my $input_gff_fn2ps      = 0;
#~ my $pout        = "-";

## my $standard_run = 1;    # get rid of !$reduced

## Get arguments (command line)
GetOptions(
    'species:s'       => \$species,
    'gff:s'           => \$input_gff_fn,
    'fasta:s'         => \$input_fas_fn,
    'sout|statsout:s' => \$sout,
    'fix_fasta'       => \$fix_fasta

      #'branch'          => \$branch_pred_flag,
      #'reduced|red'     => \$reduced,

      #'path|binpath:s'  => \$path,
      #'interactive|i' => \$interactive_flag,
      #'tenfold'       => \$tenfold,
      #'gff2ps'        => \$input_gff_fn2ps
);
my $usage =
  "Usage: $PROGRAM_NAME -species H.sapiens -gff gff_name -fasta fasta_name -sout statsfile_out -fix_fasta 1";

#~ -branch -reduced -path <executables_path>\n";

if (!($species && $input_gff_fn && $input_fas_fn && $sout && $fix_fasta))
{
    print {*STDERR} $usage and exit;
}

### Got...
### $species
### $input_gff_fn
### $input_fas_fn
### $sout
### $fix_fasta

### [<now>] Starting script at <file>[<line>]...

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
Readonly::Scalar my $train_sites_cutoff        => 1_400;
Readonly::Scalar my $train_sites_cutoff_alt    => 1_200; # changed in some part?
Readonly::Scalar my $train_sites_markov_cutoff => 5_500;
Readonly::Scalar my $backgrnd_kmer_size        => 60;
Readonly::Scalar my $backgrnd_kmer_num => 100_000;

#~ ( $totalcodingbases > 400000 && $totalnoncodingbases > 100000 )
Readonly::Scalar my $coding_bp_limit_A => 400_000;
Readonly::Scalar my $coding_bp_limit_B => 375_000;

#~ || ( $totalcodingbases > 375000 && $totalnoncodingbases > 150000 )
Readonly::Scalar my $non_coding_bp_limit_A => 100_000;
Readonly::Scalar my $non_coding_bp_limit_B => 150_000;
Readonly::Scalar my $non_coding_bp_limit_C => 35_000;

Readonly::Scalar my $multi_total_noncodingbases => 25;

Readonly::Scalar my $coding_param_X => 0.25;
Readonly::Scalar my $intron_param_X => 10;

## need to explain or just incorporate Readonly::Scalar
## BUG ??? is it ever used?
##  ## Q_FRANCISCO: condition below?
#~ || (   $totalnoncodingbases > 35000
#~ && $totalcodingbases > ( 25 * $totalnoncodingbases ) )

Readonly::Hash my %profile_params => (
    ## EXON WEIGHT PARAMETER
    ExWeightParam_ini   => -4.5,
    ExWeightParam_step  => 0.5,
    ExWeightParam_limit => -2.5,
    ## EXON/OLIGO FACTOR PARAMETER
    OligoWeight_ini   => 0.25,
    OligoWeight_step  => 0.05,
    OligoWeight_limit => 0.50,
    ## Minimum Branch Profile Distance
    ## unused 20161201
    ##    i_BranchProfDist => 7,
    ##    d_BranchProfDist => 2,
    ##    f_BranchProfDist => 9,
    ## ACCEPTOR CONTEXT
    ##    iAccCtx => 40,
    ##    dAccCtx => 10,
    ##    fAccCtx => 70,
);

Readonly::Hash my %canonical_const => (
    Donor    => [31, 'GT'],
    Acceptor => [28, 'AG'],
    Start    => [30, 'ATG'],
    ## stop 2do
    ## branch 2do
                                      );
Readonly::Hash my %my_info_thresholds => (
    donor    => 0.15,
    acceptor => 0.04,
    ATGx     => 0.15,
    ##     branch   => 0.30,
    ## stop 2do
                                         );

## end set CONSTANTS

## changing to /tmp for faster exec on clusters
## TODO: random string in dir name to avoid conflicts with other ppl runnig the script
my $work_dir = "/tmp/workdir_00_gtrain/";

my $tmp_dir      = "$work_dir/temp_00/";
my $stats_dir    = "$work_dir/stats/";
my $sites_dir    = "$work_dir/sites/";
my $plots_dir    = "$work_dir/plots/";
my $introns_dir  = "$work_dir/introns/";
my $fastas_dir   = "$work_dir/fastas/";
my $backgrnd_dir = "$work_dir/backgrnd/";
my $cds_dir      = "$work_dir/cds/";
my $geneid_dir   = "$work_dir/geneid/";
my $results_dir  = "$work_dir/results/";

my @data_dirs = (
                 $work_dir,  $tmp_dir,    $fastas_dir,  $stats_dir,
                 $sites_dir, $plots_dir,  $introns_dir, $backgrnd_dir,
                 $cds_dir,   $geneid_dir, $results_dir,
                );

create_data_dirs(@data_dirs);

## PROGRAM SPECIFIC VARIABLES (unordered...)

## my $FH_SOUT;
## my $last_bench_time;    #for benchmarking parts of the script
my $input_nodots_gff = "$work_dir/input_gff_no_dots.gff";

#~ say "\nXXX $input_nodots_gff \n";

my $genome_all_contigs_tbl = "";
my $backgrnd_kmers_fn      = "";
## my $tblseq      = "";

my @evaluation = ();
## flow control / option variables
#~ my $run_jacknife_flag = 0;
#~ my $use_allseqs_flag = 0;

#~ my $use_branch_flag   = 0;

#~ my $run_optimizeX_flag  = 0;    ## UNUSED
# $my $contig_opt_flag = 0;
## my $ext_flg          = 0;

my $train_set_gff               = "";
my $eval_set_gff                = "";
my $contigs_all_transcr_2cols   = "";
my $train_contigs_transcr_2cols = "";
my $eval_contigs_transcr_2cols  = "";
my $train_transcr_used_int      = 0;

## Golden Path stuff
my $gp_evalcontig_fa      = "";
my $gp_evalcontig_gff     = "";
my $gp_evalcontig_len_int = 0;
my $gp_evalcontig_tbl     = "";

#my $gp_eval_fa            = "";
#my $gp_eval_gff           = "";
#my $gp_eval_len_int       = 0;
#my $gp_eval_tbl           = "";

my $gp_traincontig_fa      = "";
my $gp_traincontig_gff     = "";
my $gp_traincontig_len_int = 0;
my $gp_traincontig_tbl     = "";

#my $gp_train_fa            = "";
#my $gp_train_gff           = "";
#my $gp_train_len_int       = 0;
#my $gp_train_tbl           = "";

## Weights for profiles??
my $best_ExWeightParam = 0;
my $best_OlWeight      = 0;
my $best_Acc           = 0;
my $best_Min           = 0;

## FIX used just inside a function
# my $acceptors_subprofile = "";
# my $donors_subprofile    = "";
# my $ATGx_subprofile      = "";

## my $X_geneid_sorted_gff = "";
my $seqs_eval_gff         = "";
my $train_inframestop_int = 0;
my $eval_inframestop_int  = 0;
my $locus_id              = "";
my $locus_id_new          = "";
my $intron_long_int       = 0;
my $intron_short_int      = 0;
my $intergenic_max        = 0;
my $intergenic_min        = 0;

my $eval_cds_filtered_tbl       = "";
my $eval_filtered_gff           = "";
my $eval_introns_filtered_tbl   = "";
my $eval_locusid_filtered_2cols = "";

my $train_acceptor_tbl           = "";
my $train_ATGx_tbl               = "";
my $train_cds_filtered_tbl       = "";
my $train_donor_tbl              = "";
my $train_filtered_gff           = "";
my $train_introns_filtered_tbl   = "";
my $train_locusid_filtered_2cols = "";

## my $seqs_4evaluation_listX = "";
## my $seqs_4training_listX   = "";

## my $ATGx_tbl                   = "";
## my $tmp_X_geneid_sorted_gff    = "";
my $train_2cols_seq_locusid_fn = "";
my $eval_2cols_seq_locusid_fn  = "";
my $train_transcr_lst_fn       = "";

my $train_noncanon_accept_int = 0;
my $train_noncanon_donors_int = 0;
my $train_noncanon_ATGx_int   = 0;

my $transcr_all_number = 0;
my $train_transcr_num  = 0;

#
## INITIAL CHECKS
#
## TODO  fasta / gff file accessible?
#~ sub validate_input_fasta
#~ {
#~ return 1;
#~ }
#~ sub validate_input_gff
#~ {
#~ return 1;
#~ }
#~ sub select_test_eval
#~ {
#~ return 1;
#~ }
#~ sub extract_CDS
#~ {
#~ return 1;
#~ }
#~ sub extract_introns
#~ {
#~ return 1;
#~ }

#~ sub extract_donors
#~ {
#~ return 1;
#~ }
#~ sub extract_acceptors
#~ {
#~ return 1;
#~ }
#~ sub extract_ATGx
#~ {
#~ return 1;
#~ }
## TODO 2. limits:
## 2a. >= 500 genes in gff
## PYTHON INLINE START

#~ use Inline Python => <<END;

#~ def Foo():
#~ class Bar:
#~ def __init__(self):
#~ print( "new Bar()")
#~ def tank(self):
#~ return 10
#~ return Bar()

#~ END

#~ use Inline Python => <<'ENDPY';

#~ from pygeneid import check_fasta

#~ def JAxH(x):
#~ return "Just Another %s Hacker" % x

#~ ENDPY

#~ print JAxH('Inline'), "\n";

#~ my $o = Foo();
#~ print $o->tank(), "\n";
#~ die;
#~ use Inline Python => << 'PYEND';
#~ from pygeneid import check_fasta

#~ ##print("from python1")
#~ PYEND

#~ my $headers_fasta_seq;
#~ $headers_fasta_seq  = check_fasta($input_fas_fn);

#~ print {*STDERR} "\n PYTHON The user has provided $headers_fasta_seq  genomic sequences\n";

#~ die;
#~ ## PYTHON INLINE END
#
## Common tasks
#
## sanity check

## CREATE BLANK PARAMETER FILE##
### [<now>] template parameter file at <file>[<line>]
my $param        = Geneid::Param->new($species);
my $new_param_fn = "$work_dir/$species.geneid.param";

#set isochores to 1
$param->numIsocores(1);
$param->isocores([Geneid::Isocore->new()]);

## END CREATING A PARAMETER FILE REGARDLESS OF WHETHER THE TRAINING IS COMPLETE OR REDUCED

normal_run();

sub normal_run
{
    ## BUG, simplification needed
    ## simplification
    #    my $optymize_type = "";
    #    $optymize_type       = "contig";
    #    $contig_opt_flag = 1;

    $genome_all_contigs_tbl = $work_dir . $species . ".genomic_all_contigs.tbl";
    my $my_command =
      "./bin/fas_to_tbl.py $input_fas_fn $genome_all_contigs_tbl $fastas_dir";

    # $fastas_dir";
    run($my_command);

    #fasta_2_tbl( $input_fas_fn, $genome_all_contigs_tbl );

    run("sort -o $genome_all_contigs_tbl $genome_all_contigs_tbl");

    #tbl_2_single_fastas( $genome_all_contigs_tbl, $fastas_dir);
    #write_sizes_from_tbl_fn($genome_all_contigs_tbl);
    $my_command = "cp $input_gff_fn $input_nodots_gff";
    say "$my_command";
    run($my_command);

    #~ my $old_option = 0;

    #~ ## Convert fasta to tabular format
    #~ ## Fasta process to sub later
    #~ ## my $t0 = Benchmark->new;

    #~ if ($old_option)
    #~ {
    #~ print "not yet\n";
    #~ }

    #     print {*STDERR}
    #       "\nConverting genomics fasta file ($input_fas_fn) to tabular format\n";

    #     my $genomic_temp_tbl = $work_dir . $species . ".genomic.tmp.tbl";
    #     fasta_2_tbl( $input_fas_fn, $genomic_temp_tbl );
    #     run("sort -o $genomic_temp_tbl $genomic_temp_tbl");

    #     print {*STDERR} "actg to ACTG conversion of input fasta \n";
    #     my $tblcaps = "";

    #     open(
    #         my $FH_LOCID,
    #         "-|",
    # "gawk '{gsub(/_/,\"\",\$1);gsub(/\\./,\"\",\$1);print \$1, toupper(\$2)}' $genomic_temp_tbl "
    #     ) or croak "Failed here";
    #     while (<$FH_LOCID>) {
    #         $tblcaps .= $_;
    #     }
    #     close $FH_LOCID;

    #     chomp $tblcaps;
    #     $genome_all_contigs_tbl = $work_dir . $species . ".genomic.tbl";
    #     open( my $FH_FOUT_caps, ">", "$genome_all_contigs_tbl" )
    #       or croak "Failed here";
    #     print {$FH_FOUT_caps} "$tblcaps";
    #     close $FH_FOUT_caps;

    # ## place genomic sequences in "fastas_$species" directory
    #     print {*STDERR} "move genomic sequences into \"$fastas_dir\" directory\n";
    #     print {*STDERR} "(also transfer genomic fasta length info)\n\n";
    # ## do not create fastas in diretory if they are already created and their number corresponds to the number of sequences in thr array
    # ## CONVERT GENOMICS FASTA TO MULTI FASTA AND PLACE THEM IN APPROPRIATE DIRECTORY

    #     print {*STDERR}
    # "Convert $genome_all_contigs_tbl to multiple genomic fastas and place them in $fastas_dir:\n";

    #     tbl_2_single_fastas( $genome_all_contigs_tbl, $fastas_dir );
    #     print {*STDERR}
    # "\n\nConversion of $genome_all_contigs_tbl to multiple genomic fastas completed..\n\nAdd fasta sequence length information to same directory\n\n";
    #     write_sizes_from_tbl_fn($genome_all_contigs_tbl);

    # #

    ## get locus_id file only first time pipeline is run for a given species #ALL GENE MODELS
    ##  ## Q_FRANCISCO: can we assume some sane GFF/GFT format as an input?  Why _ and "."?

    #     if ($old_option) {
    #         print {*STDERR}
    #           "\nEliminate undesirable (_ and .) characters from $input_gff_fn\n";

    #         my $filtergff = "";

    #         $my_command =
    # "gawk '{OFS=\"\\t\"}{gsub(/\\./,\"\",\$1);gsub(/\\./,\"\",\$9);gsub(/_/,\"\",\$0);print}' $input_gff_fn";
    #         $filtergff = capture($my_command);

    #         open( my $FH_FOUT, ">", "$input_nodots_gff" ) or croak "Failed here";
    #         print {$FH_FOUT} "$filtergff";
    #         close $FH_FOUT;

    #     }
    #     else {
    #         $input_nodots_gff = $input_gff_fn;
    #     }
    ### Obtain list of genomic sequences at <file>[<line>]...

    $my_command = "gawk '{print \$1,\$9}' $input_nodots_gff | sort | uniq ";
    $locus_id   = capture($my_command);

    # my $FH_FOUT;
    say "\n TTT got here TTT\n";
    $contigs_all_transcr_2cols = $work_dir . $species . "_locus_id";
    open(my $FH_FOUT, ">", "$contigs_all_transcr_2cols") or croak "Failed here";
    print {$FH_FOUT} "$locus_id";
    close $FH_FOUT;

    ## number of gene models TOTAL
    $my_command =
      " gawk '{print \$2}' $contigs_all_transcr_2cols | sort | uniq | wc -l";
    $transcr_all_number = capture($my_command);

    chomp $transcr_all_number;
    ## number of genomic sequences TOTAL
    $my_command =
      "gawk '{print \$1}' $contigs_all_transcr_2cols | sort | uniq | wc -l";
    my $total_genomic = capture($my_command);

    chomp $total_genomic;

    print {*STDERR}
      "\nThe gff file ($input_nodots_gff) contains a total of $total_genomic genomic sequences and $transcr_all_number gene models\n";

    ## get a list of genes TOTAL
    ### Obtain list of all transcripts at <file>[<line>]...

    $my_command = "gawk '{print \$9}' $input_nodots_gff | sort | uniq ";
    my $transcr_list_tmp = capture($my_command);

    my $transcr_all_list_fn = $work_dir . $species . "_all_seqs.lst";
    open($FH_FOUT, ">", "$transcr_all_list_fn") or croak "Failed here";
    print {$FH_FOUT} "$transcr_list_tmp";
    close $FH_FOUT;

    if ($transcr_all_number >= $train_loci_cutoff)
    {

        $train_transcr_num = int($train_fraction * $transcr_all_number);

        print {*STDERR}
          "\nA subset of $train_transcr_num sequences (randomly chosen from the $transcr_all_number gene models) was used for training\n";
        ## DEBUG KEEP !!! shuf => random select
        ## head -$train_transcr_num just the first ones
        #my $my_command =           "shuf --head-count=$train_transcr_num $contigs_all_transcr_2cols | sort ";

        $my_command =
          "head --lines=$train_transcr_num $contigs_all_transcr_2cols ";

        $locus_id_new = capture($my_command);

        $train_contigs_transcr_2cols =
          $work_dir . $species . "_train_setaside80.2cols";
        open($FH_FOUT, ">", "$train_contigs_transcr_2cols")
          or croak "Failed here";
        print {$FH_FOUT} "$locus_id_new";
        close $FH_FOUT;

        ## ASSUMING USER SELECTED TO SET ASIDE SEQUENCES FOR EVALUATION (20%)
        ## 2016.12.10 $my_command = "gawk '{print \$2}' $train_contigs_transcr_2cols | wc -l";
        ## 2016.12.10  $train_transcr_used_num = capture($my_command);
        ## 2016.12.10 chomp $train_transcr_used_num;
        $train_transcr_used_int = $train_transcr_num;

        #~ my $t1 = Benchmark->new;
        #~ my $td = timediff( $t1, $t0 );
        #~ print "\nTTT the code took t0->t1:", timestr($td), "\n";
        #~ $last_bench_time = $t1;

        #
        ## gff for training subset
        ##
        my $gff_4_training = "";

        print {*STDERR}
          "\nThe new training gff file includes $train_transcr_used_int gene models (80% of total seqs)\n";
        ## ??? BUG ???
        $my_command =
          "gawk '{print \$2\"\$\"}' $train_contigs_transcr_2cols | sort | uniq | egrep -wf - $input_nodots_gff";
        $gff_4_training = capture($my_command);

        $train_set_gff = $work_dir . $species . "_train_setaside80.gff";
        open($FH_FOUT, ">", "$train_set_gff") or croak "Failed here";
        print {$FH_FOUT} "$gff_4_training";
        close $FH_FOUT;

        print {*STDERR} "\nObtain list of training genes\n\n";

        #my $train_transcr_list_tmp = "";

        $my_command = "gawk '{print \$9}' $train_set_gff | sort | uniq ";
        my $train_transcr_list_tmp = capture($my_command);

        $train_transcr_lst_fn = $work_dir . $species . "train_setaside80.lst";
        open($FH_FOUT, ">", "$train_transcr_lst_fn") or croak "Failed here";
        print {$FH_FOUT} "$train_transcr_list_tmp";
        close $FH_FOUT;

        #
        ## new locus_id for evaluation test set
        #
        my $locus_id_eval = "";

        $my_command =

          #"gawk '{print \$0\"\$\"}' $train_transcr_lst_fn | egrep -vwf - $contigs_all_transcr_2cols";
          "grep -vwf $train_transcr_lst_fn $contigs_all_transcr_2cols";

        $locus_id_eval = capture($my_command);
        chomp $locus_id_eval;

        $eval_contigs_transcr_2cols =
          $work_dir . $species . "_evaluation_setaside20.2cols";
        open($FH_FOUT, ">", "$eval_contigs_transcr_2cols")
          or croak "Failed here";
        print {$FH_FOUT} "$locus_id_eval";
        close $FH_FOUT;

        ##
        ## gff for evaluation test set
        #

        $my_command =
          "gawk '{print \$2\"\$\"}' $eval_contigs_transcr_2cols | sort | uniq | egrep -wf - $input_nodots_gff | gawk '{ print \$9}' | sort | uniq | wc -l";
        ## ??? BUG this is a number...
        $seqs_eval_gff = capture($my_command);

        #chomp $input_gff_fnseqseval;

        print {*STDERR}
          "The evaluation gff file includes $seqs_eval_gff gene models (20% of total seqs)\n\n";

        $my_command =
          "gawk '{print \$2\"\$\"}' $eval_contigs_transcr_2cols | sort | uniq | egrep -wf - $input_nodots_gff ";
        my $gff_4_evaluation = capture($my_command);

        $eval_set_gff = $work_dir . $species . "_evaluation_setaside20.gff";
        open($FH_FOUT, ">", "$eval_set_gff") or croak "Failed here";
        print {$FH_FOUT} "$gff_4_evaluation";
        close $FH_FOUT;

    }    # seqs > 500
    ##LOOP IF WE HAVE FEWER THAN 500 SEQUENCES

    else
    {    # seqs < $train_loci_cutoff
        ## BUG we do not do jacknife anyway here
        croak "we do not have >= $train_loci_cutoff sequences, quitting now";
    }    # seqs < 500

    #~ if (!$use_allseqs_flag)
    #~ {    ##SET SEQS FOR EVAL AND TRAINING (SUBSETS)

    ### Not using all genes for training <file>[<line>]...

    #~ ## not used
    #~ #  "\nConvert general gff2 to geneid-gff format  NOT_USE_ALL_SEQS \n\n";
    #~ # XXX function name
    #~ }
    ## Convert general gff2 to geneid gff format
    ## extract and check cds and intron sequences. Remove inframe stops and check all seqs start with ATG and end with STOP
    ## TRAIN
    #~ my $t1 = Benchmark->new;
    #~ my $td = timediff($t1, $t0);
    #~ print "\nTTT the code took t0->t1:",timestr($td),"\n";
    #~ $last_bench_time = $t1;

    ## XXXXXX extracted lines start
    #
    ## CREATE FASTAS CDS; INTRON, SITES DIRs WITHIN PATH (ONLY FIRST TIME)
    {
        ## BUG => there is no need to convert gff to geneid format each time we run

        $train_2cols_seq_locusid_fn =
          gff_2_geneidgff_mock($train_set_gff, $species, ".train");
        ### [<now>] at <file>[<line>]...
        ### $train_2cols_seq_locusid_fn

        (
         $train_cds_filtered_tbl,       $train_introns_filtered_tbl,
         $train_locusid_filtered_2cols, $train_filtered_gff,
         $train_inframestop_int
        )
          = @{
            extract_cds_introns(
                                $train_2cols_seq_locusid_fn,
                                $train_contigs_transcr_2cols,
                                ".train"
                               )
             };

        #  print {*STDERR} " OUTSIDE EXTRACTCDSINTRON outgff: $train_filtered_gff\noutlocus_id: $out_locus_id_X\n";
        ## TRAIN
        ## EVAL
        $eval_2cols_seq_locusid_fn =
          gff_2_geneidgff_mock($eval_set_gff, $species, ".eval");

        ### [<now>] at <file>[<line>]...
        ### $eval_2cols_seq_locusid_fn
        ### $eval_contigs_transcr_2cols\n";
        ## 2016.12.10 not used??
        (
         $eval_cds_filtered_tbl, $eval_introns_filtered_tbl,    #-- not used
         $eval_locusid_filtered_2cols,                          #-- not used
         $eval_filtered_gff, $eval_inframestop_int
        )
          = @{
            extract_cds_introns(
                                $eval_2cols_seq_locusid_fn,
                                $eval_contigs_transcr_2cols,
                                ".eval"
                               )
             };

        #EVAL

    }

    #
    ### [<now>] Get background kmers next <file>[<line>]...
    ### $backgrnd_kmer_num
    ### $backgrnd_kmer_size

    $backgrnd_kmers_fn = $work_dir . $species . "_background.tbl";
    get_background_kmers(
                         $backgrnd_kmer_size,     $input_fas_fn,
                         $genome_all_contigs_tbl, $backgrnd_kmer_num,
                         $backgrnd_kmers_fn
                        );
    ## extract and check splice sites and start codon. Use only canonical info #IN SEQUENCES USED IN TRAINING
    ### [<now>] at <file>[<line>]...
    ### $train_filtered_gff
    ### $train_locusid_filtered_2cols
    (
     $train_donor_tbl,    $train_noncanon_donors_int,
     $train_acceptor_tbl, $train_noncanon_accept_int,
     $train_ATGx_tbl,     $train_noncanon_ATGx_int
    )
      = @{extract_sites($train_filtered_gff, $train_locusid_filtered_2cols)};

    ## prepare sequences for optimization of newly developed parameter file (TRAIN)

    ### [<now>] Convert gff to gp (golden-path-like)format <file>[<line>]...
    ### artificial contig - concatenated sequences - approx. 800bp linkers...

    (
     $gp_traincontig_gff, $gp_traincontig_fa,
     $gp_traincontig_tbl, $gp_traincontig_len_int
    )
      = @{process_seqs_4opty($train_filtered_gff, ".train", 1)};

    #~ print {*STDERR}
    #~ "\nConvert gff to gp (golden-path-like)format (training set for later optimization -400-nt flanked sequences)\n";
    #~ ($gp_train_gff, $gp_train_fa, $gp_train_tbl, $gp_train_len_int) =
    #~ @{process_seqs_4opty($train_filtered_gff, ".train", 0)};
    #~ print {*STDERR} "$gp_train_gff";

    ## prepare test set for evaluation of newly developed parameter file (EVAL)

    #~ print {*STDERR}
    #~ "\nConvert gff to gp (golden-path-like)format (400-nt flanking)(test set for evaluation of new parameter file)\n";
    #~ ($gp_eval_gff, $gp_eval_fa, $gp_eval_tbl, $gp_eval_len_int) =
    #~ @{process_seqs_4opty($eval_filtered_gff, ".eval", 0)};
    #~ print {*STDERR} "DONE\n";

    ### [<now>] at <file>[<line>]...

    (
     $gp_evalcontig_gff, $gp_evalcontig_fa,
     $gp_evalcontig_tbl, $gp_evalcontig_len_int
    )
      = @{process_seqs_4opty($eval_filtered_gff, ".eval", 1)};

    #~ my $t3 = Benchmark->new;
    #~ my $td = timediff( $t3, $last_bench_time );
    #~ print "\nTTT the code took t2->t3:", timestr($td), "\n";
    #~ $last_bench_time = $t3;

    # my (
    #     $donor_start,  $donor_end,  $acceptor_start,
    #     $acceptor_end, $ATGx_start, $ATGx_end
    # ) = compute_sites_pictogram();

    my ($donor_start, $donor_end) =
      compute_matrices_4sites($train_donor_tbl, 'donor');
    my ($acceptor_start, $acceptor_end) =
      compute_matrices_4sites($train_acceptor_tbl, 'acceptor');
    my ($ATGx_start, $ATGx_end) =
      compute_matrices_4sites($train_ATGx_tbl, 'ATGx');

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
       )
      = @{
        derive_coding_potential(
                                $train_cds_filtered_tbl,
                                $train_introns_filtered_tbl
                               )
         };

    #add markov matrices to the parameter file
    if (!defined @{$param->isocores}[0]->Markov_order($markov_model))
    {
        croak "error in setting Markov_order\n";
    }
    if (!defined @{$param->isocores}[0]
        ->Markov_Initial_probability_matrix($markov_mod_ini))
    {
        croak "error in setting Markov_Initial_probability_matrix\n";
    }
    if (!defined @{$param->isocores}[0]
        ->Markov_Transition_probability_matrix($markov_mod_trans))
    {
        croak "error in setting Markov_Transition_probability_matrix\n";
    }
    ##

    ## BUG do not remove, misleading name of a function: WriteStatsFile
    ($intron_short_int, $intron_long_int, $intergenic_min, $intergenic_max) =
      calc_stats(
        ## $species,                   $sout,

        $train_introns_filtered_tbl, $train_cds_filtered_tbl,
        $train_filtered_gff,         $train_inframestop_int,
        $eval_inframestop_int,       $train_transcr_used_int,
        $train_noncanon_donors_int,  $train_noncanon_accept_int,
        $train_noncanon_ATGx_int,    $markov_model,
        $total_coding,               $total_noncoding,
        $donor_start,                $donor_end,
        $acceptor_start,             $acceptor_end,
        $ATGx_start,                 $ATGx_end,

        #~ 0,
        #~ 0,
        #~ 0,
        #~ #$use_allseqs_flag
      );

    ### [<now>] Intron/intergenic clculations <file>[<line>]...
    ### $intron_short_int
    ### $intron_long_int
    ### $intergenic_min
    ### $intergenic_max

    ##
    ### [<now>] Write preliminary non-optimized parameter file at <file>[<line>]...

    $new_param_fn = "$work_dir/$species.geneid.param";
    $param->writeParam($new_param_fn);

    # if reduced training (non-default) do not calculate any of the above ALL OF THE ABOVE MUST BE RUN ONLY FIRST TIME GENEID IS TRAINED FOR A GIVEN SPECIES
    #EVERYTHING BELOW WILL BE RUN EVERYTIME THE TRAINING PIPELINE IS RUN WHETHER "REDUCED" OR "FULL"

    ### [<now>] Optimizing new parameter file...

    #~ $run_jacknife_flag   = 0;

    ## BUG settings need to be set forward. Also these are numbers.

    ## OPTIMIZATION FUNCTION NO BRANCH
    my $array_ref         = "";
    my $OligoWeight_ini   = $profile_params{'OligoWeight_ini'};
    my $ExWeightParam_ini = $profile_params{'ExWeightParam_ini'};

    ## OPTIMIZATION FUNCTIONS
    #
    #    if (!$contig_opt_flag)
    #    {
    #        print {*STDERR} "\n DEBUG: NOT CONTIG OPT\n";
    #        croak "bad bug!\n";
    #    }    #end if
    #~ if ( !$contig_opt_flag ) {
    #~ print {*STDERR} "\n DEBUG: NOT CONTIG OPT\n";
    #~ @evaluation = @{
    #~ parameter_optimize( $gp_train_fa, $gp_train_gff, $new_param_fn, 0,
    #~ 0, 0, 0,
    #~ $IeWF, $deWF, $FeWF, $IoWF, $doWF, $FoWF, 0, 0, 0, 0, 0, 0 )
    #~ };

    #~ ( $best_IeWF, $best_IoWF, $best_Acc, $best_Min, $array_ref ) = @{
    #~ BuildOptimizedParameterFile( \@evaluation, $use_branch_flag, 0, 0,
    #~ 0 )
    #~ };

    #~ }
    ## elsif ($contig_opt_flag)

    ### [<now>] contig optimisation...
    ## too many parameters in old code
    #        @evaluation = @{
    #            parameter_optimize( $gp_traincontig_fa, $gp_traincontig_gff,
    #                $new_param_fn, 0, 0, 0, 0, $IeWF, $deWF, $FeWF, $IoWF, $doWF,
    #                $FoWF, 0, 0, 0, 0, 0, 0 )
    #        };

    @evaluation = @{
        parameter_optimize(
                           $gp_traincontig_fa,
                           $gp_traincontig_gff,
                           $new_param_fn,
                           %profile_params
                          )
                   };

    ($best_ExWeightParam, $best_OlWeight, $best_Acc, $best_Min, $array_ref) = @{
        get_opt_paramfile(\@evaluation)
          ##, $use_branch_flag, 0, 0, 0 )
    };

    my @evaluation_init = @{$array_ref};
    my @evaluation_test = ();

    ##
    ## EVALUATE PERFORMANCE OF NEW PARAMETER FILE ON TEST SET (IF PRESENT)
    ##

    my $param_opt_fn = "$species.geneid.optimized.param";

    #~ if (!$use_allseqs_flag)
    #~     {
    ##my $FH_SOUT;
    open(my $FH_SOUT, ">", "$work_dir/$species.use_NOT_allseqs.log");

    #print {*STDERR} "CHECK EVALUATE: $gp_eval_fa, $gp_eval_gff, $param_opt_fn\n";
    #
    #    if (!$contig_opt_flag)
    #    {
    #
    #        @evaluation_test = @{
    #            parameter_evaluate(
    #                $gp_eval_fa,      $gp_eval_gff, $param_opt_fn,
    #                $OligoWeight_ini, $ExWeightParam_ini
    #            )
    #            };
    #
    #    }    #end if !$contig_opt_flag

    ##elsif ($contig_opt_flag)
    ##{

    @evaluation_test = @{
        parameter_evaluate(
                           $gp_evalcontig_fa, $gp_evalcontig_gff, $param_opt_fn,
                           $OligoWeight_ini,  $ExWeightParam_ini
                          )
                        };
    ##}

    ##~ if (!$use_branch_flag)
    #~ {

    print {*STDERR}
      "\nPerformance of new optimized parameter file on test set:\n\n"
      . join("\t", @evaluation_init[2 .. $#evaluation_init]), "\n";

    #~ }
    #~ elsif ($use_branch_flag)
    #~ {

    #~ print {*STDERR}
    #~ "\nPerformance of new optimized parameter file on test set:\n\n"
    #~ . join("\t", @evaluation_init[4 .. $#evaluation_init]), "\n";

    #~ }

    print {*STDERR} join("\t", @evaluation_test), "\n\n";

    print {$FH_SOUT} join("\t", @evaluation_test), "\n\n";
    close $FH_SOUT;

    #~ }    # if NOT using all seqs for training

    return 1;
} ## end normal run

#
## END OF MAIN PORTION OF SCRIPT
#

sub extract_cds_introns ($my_nodots_gff, $my_contig_transcr_2cols, $type)
{

    ### [<now>] Running extract_cds_introns at <file>[<line>]...
    ### Got: $type
    ### $my_contig_transcr_2cols
    ### $my_nodots_gff

    my $out_cds_intron_gff =
      $work_dir . $species . "_" . $type . ".cds_intron_gff";
    open(my $FH_LOCUS, "<", "$my_contig_transcr_2cols")
      or croak "Failed here";
    print {*STDERR} "$my_contig_transcr_2cols and $my_nodots_gff\n";
    my $count = 0;
    while (<$FH_LOCUS>)
    {
        my ($id_genomic, $id_gene) = split;
        my $tmp_single_gene_gff = "$tmp_dir/$id_gene.gff";
        run(" egrep -w '$id_gene\$' $my_nodots_gff | sort -k4,5n > $tmp_single_gene_gff"
           );
        ## NEW code simplify
        run(" egrep -w '$id_gene\$' $my_nodots_gff | sort -k4,5n >> $out_cds_intron_gff"
           );
        ## POTENTIAL BUG, split commands below

        #~ my $FH_ssgff_A    = File::Temp->new();
        #~ my $fname_ssgff_A = $FH_ssgff_A->filename;
        ## CDS
        my $fname_ssgff_A = "$cds_dir/${species}_ssgff_A.tmp.fa";
        my $my_command =
          "./bin/ssgff -cE $fastas_dir/$id_genomic $tmp_single_gene_gff >  $fname_ssgff_A";
        run($my_command);
        $my_command =
          "cat $fname_ssgff_A | sed -e 's/:/_/' -e 's/ CDS//' >> $cds_dir/${species}${type}.cds.fa ";
        run($my_command);

        #` ./bin/ssgff -cE $work_dir/fastas_$species/$id_genomic $tmp_dir/$id_gene.gff | sed -e 's/:/_/' -e 's/ CDS//' >> $work_dir/cds/${species}${type}.cds.fa `;

        ## INTRONS
        #~ my $FH_ssgff_B    = File::Temp->new();
        #~ my $fname_ssgff_B = $FH_ssgff_B->filename;
        my $fname_ssgff_B = "$introns_dir/${species}_ssgff_B.tmp.fa";
        $my_command =
          "./bin/ssgff -iE $fastas_dir/$id_genomic $tmp_single_gene_gff > $fname_ssgff_B";

        #say "\n$my_command\n";
        run($my_command);
        $my_command =
          "cat $fname_ssgff_B | sed -e 's/:/_/' -e 's/ Intron.*//' >> $introns_dir/${species}${type}.intron.fa";

        #say "\n$my_command\n";
        run($my_command);

        #` ./bin/ssgff -iE $work_dir/fastas_$species/$id_genomic $tmp_dir/$id_gene.gff | sed -e 's/:/_/' -e 's/ Intron.*//' >> $work_dir/intron/${species}${type}.intron.fa `;
        $count++;
        print {*STDERR} "$count ..";
    }
    close $FH_LOCUS;
### [<now>] Finished extracting introns/CDS...
    #    print {*STDERR} "DONE\n";

    # #tabulate CDS and INTRON SEQUENCES

    ### [<now>] Create tabular format of CDS and intron sequences at <file>[<line>]...
    ### $type

    #    print {*STDERR}
    #        "\nCreate tabular format of CDS and INTRON sequences for $type sequences\n";

    ## CDS
    my $cds_tmp_fa = $cds_dir . ${species} . "$type" . ".cds.fa";

    ### $cds_tmp_fa
    #    print {*STDERR} "$cds_tmp_fa\n\n";

    ### [<now>] cds fasta_2_tbl at <file>[<line>]...
    my $cds_tmp_tbl = $cds_dir . ${species} . "$type" . ".cds.tbl";
    fasta_2_tbl($cds_tmp_fa, $cds_tmp_tbl);

    #    print {*STDERR} "cds tabular file created for $type sequences \n";

    # ##INTRON
    my $intron_tmp_fa  = $introns_dir . ${species} . "$type" . ".intron.fa";
    my $intron_tmp_tbl = $introns_dir . ${species} . "$type" . ".intron.tbl";
    ### [<now>] intron fasta_2_tbl at <file>[<line>]...
    fasta_2_tbl($intron_tmp_fa, $intron_tmp_tbl);

    # INTRONS LARGER THAN 0 ONLY
    #
    #    my $intron_nonzero_tmp = "";
    #    my $my_command =
    #      "gawk '{if(length(\$2)>0){print \$1,\$2}}' $intron_tmp_tbl ";
    #    $intron_nonzero_tmp = capture($my_command);
    #
    #    my $intron_nonzero_tbl =
    #      $work_dir . $species . "$type" . ".intron_positivelength.tbl";
    #
    #    open( my $FH_FOUT, ">", "$intron_nonzero_tbl" ) or croak "Failed here";
    #    print {$FH_FOUT} "$intron_nonzero_tmp";
    #    close $FH_FOUT;
    #
    #    print {*STDERR}
    #      "intron tabular file created with introns with more than 0 nucleotides\n";
    # GET LIST OF SEQUENCES WITH LENGTH >0 and EXCLUDE FROM CDS/locus_id/gff FILES SEQUENCES WITH INTRONS WITH 0 LENGTH
    #    my $intron_zero = "";
    #    $my_command =
    #"gawk '{if(length(\$2)==0){print \$1}}' $intron_tmp_tbl | sed 's/\\(.*\\)\\..*/\\1\\_/' | sort | uniq ";
    #    $intron_zero = capture($my_command);
    #
    #    my $temp_all_intron_zero_list =
    #      $work_dir . $species . "$type" . ".intron_zerolength.list";
    #
    #    open( $FH_FOUT, ">", "$temp_all_intron_zero_list" );
    #    print {$FH_FOUT} "$intron_zero";
    #    close $FH_FOUT;
    #
    #    my $intron_zero2 = "";
    #    $my_command =
    #"gawk '{if(length(\$2)==0){print \$1}}' $intron_tmp_tbl | sed 's/\\(.*\\)\\..*/\\1/' | sort | uniq ";
    #    $intron_zero2 = capture($my_command);
    #
    #    my $temp_all_intron_zero_list2 =
    #      $work_dir . $species . "$type" . ".intron_zerolength.list2";
    #
    #    open( $FH_FOUT, ">", "$temp_all_intron_zero_list2" );
    #    print {$FH_FOUT} "$intron_zero2";
    #    close $FH_FOUT;
    #
    # FILTER SEQUENCES WITH 0 SIZE INTRONS FROM CDS!
    #
    #    $my_command = "egrep -vf $temp_all_intron_zero_list $cds_temp_tbl ";
    #    my $cds_nozero_tbl = capture($my_command);
    #
    #    my $cds_all_nozero_tbl = $work_dir . $species . "$type" . ".cds_nozero.tbl";
    #
    #    open( $FH_FOUT, ">", "$cds_all_nozero_tbl" );
    #    print {$FH_FOUT} "$cds_nozero_tbl";
    #    close $FH_FOUT;
    # ENSURE LOCUSID DOES NOT CONTAIN SEQUENCES WITH 0 SIZE INTRONS
    #    $my_command = "egrep -vwf $temp_all_intron_zero_list2 $locus_id ";
    #    my $locusid_nozero = capture($my_command);
    #
    #    my $temp_locus_id_nozero =
    #      $work_dir . $species . "$type" . "_locus_id_nozero";
    #    open( $FH_FOUT, ">", "$temp_locus_id_nozero" );
    #    print {$FH_FOUT} "$locusid_nozero";
    #    close $FH_FOUT;
    # ENSURE GFF DOES NOT CONTAIN SEQUENCES WITH 0 SIZE INTRONS
    #
    #    my $gffnozero = "";
    #    $my_command = "egrep -vwf $temp_all_intron_zero_list2 $my_nodots_gff ";
    #    $gffnozero  = capture($my_command);
    #
    #    my $exCI_temp_nonzero_gff =
    #      $work_dir . $species . "$type" . ".non_zero.gff";
    #
    #    open( $FH_FOUT, ">", "$exCI_temp_nonzero_gff" ) or croak "Failed here";
    #    print {$FH_FOUT} "$gffnozero";
    #    close $FH_FOUT;
    #
    #    #    rmtree([ "$path/cds/" ]);
    #    #    rmtree([ "$path/intron/" ]);
    #
    # Convert sequences to protein format and check for in-frame stops
    #    print {*STDERR}
    #"\nConvert sequences to protein format and check for in-frame stops and for proteins not starting with an M or not ending with a STOP\n\n";
    #
    # SHOWS WHERE GENETIC CODE FILE IS LOCATED AND ITS NAME
    #
    #    my $temp_all_protein = $work_dir . $species . "$type" . ".protein";
    #
    ## $temp_all_protein = translate_2_protein($genetic_code,$temp_cds,$temp_all_protein);
    #    $temp_all_protein =
    #      translate_2_protein( $genetic_code, $cds_all_nozero_tbl,
    #        $temp_all_protein );
    #
    #    $my_command =
    #"gawk '{print \$2,\$1}' $temp_all_protein | egrep '[A-Z]\\*[A-Z]\|^[^M]\|[^\\*] ' | gawk '{print \$2}' | wc -l";
    #    ## BUG reports just the number
    #    my $inframestops = capture($my_command);
    #    chomp $inframestops;
    #
    ##~ #my $inframe_Xstops =`gawk '{print \$2,\$1}' $temp_all_protein | egrep '[A-Z]\\*[A-Z]\|^[^M]\|[^\\*] ' | gawk '{print \$2}' | wc | gawk '{print \$1}'`;
    ##~ chomp $inframe_Xstops;
    #
    #    print {*STDERR}
    #"\n\nWe found $inframestops sequences with in-frame stop signals/not starting with a methionine or not ending with a canonical stop codon \n\n";
    #
    # IF INFRAME
    #    print $temp_all_protein;
    #
    #    if ($inframestops) {
    #        my $inframe = "";
    #        my @inframe = ();
    #
    #        ## XXX may not work XXX
    #        $my_command =
    #"gawk '{print \$2,\$1}' $temp_all_protein | egrep '[A-Z]\\*[A-Z]|^[^M]|[^\\*]' | gawk '{print \$2}' | sort | uniq ";
    #        @inframe = capture($my_command);
    #
    #        foreach my $line (@inframe) {
    #            my (@frame) = split "_", $line;
    #            my $first = $frame[0];
    #            $inframe .= "$first\n";
    #        }
    #
    #        my $inframe_protein =
    #          $work_dir . $species . "$type" . "_INframe_NoMethionine_NoSTOP";
    #        open( $FH_FOUT, ">", "$inframe_protein" ) or croak "Failed here";
    #        print {$FH_FOUT} "$inframe";
    #        close $FH_FOUT;
    # REMOVE SEQUENCES WITH IN-FRAME STOPS FROM ORIGINAL CDS / INTRON / LOCUS_ID /GFF FILES AND PRINT NEW FILES
    #        print {*STDERR}
    #"\nremove sequences with in-frame stop signals from cds/intron files\n\n";
    #
    #        $my_command =
    #"sed 's/\\(.*\\)/\\1_/g' $inframe_protein | egrep -vf - $cds_all_nozero_tbl";
    #        my $cdstbl2 = capture($my_command);
    #
    #        my $temp_all_cds2 = $work_dir . $species . "$type" . ".cds_filter1.tbl";
    #
    #        open( $FH_FOUT, ">", "$temp_all_cds2" );
    #        print {$FH_FOUT} "$cdstbl2";
    #        close $FH_FOUT;
    #
    #        my $introntbl2 = "";
    #
    #        $my_command =
    #"sed 's/\\(.*\\)/\\1\.i/g' $inframe_protein | egrep -vf - $intron_nonzero_tbl ";
    #        $introntbl2 = capture($my_command);
    #
    #        my $temp_all_intron2 =
    #          $work_dir . $species . "$type" . ".intron_filter1.tbl";
    #
    #        open( $FH_FOUT, ">", "$temp_all_intron2" );
    #        print {$FH_FOUT} "$introntbl2";
    #        close $FH_FOUT;
    #
    #        $my_command =
    #"sed 's/\\(.*\\)/\\1\$/g' $inframe_protein | egrep -vf - $temp_locus_id_nozero ";
    #        my $new_locus_id_filter1 = capture($my_command);
    #
    #        my $temp_locus_id_new2 =
    #          $work_dir . $species . "$type" . "_locus_id_filter_noinframe";
    #
    #        open( $FH_FOUT, ">", "$temp_locus_id_new2" );
    #        print {$FH_FOUT} "$new_locus_id_filter1";
    #        close $FH_FOUT;
    #
    #        #my $gffnew = "";
    #        $my_command =
    #"sed 's/\\(.*\\)_.*/\\1\$/g' $inframe_protein | egrep -vf - $exCI_temp_nonzero_gff ";
    #        my $gffnew = capture($my_command);
    #
    #        my $temp_newgff = $work_dir . $species . "$type" . ".noinframe.gff";
    #
    #        open( $FH_FOUT, ">", "$temp_newgff" ) or croak "Failed here";
    #        print {$FH_FOUT} "$gffnew";
    #        close $FH_FOUT;
    #
    ##
    #        return [
    #            $temp_all_cds2, $temp_all_intron2, $temp_locus_id_new2,
    #            $temp_newgff,   $inframestops
    #        ];

    ##    }
    #    else {    ## ??? END IF THERE ARE INFRAME STOPS
    my $my_cds_filtered_tbl     = $cds_tmp_tbl;
    my $my_intron_filtered_tbl  = $intron_tmp_tbl;
    my $my_locusid_filtered_tbl = $my_contig_transcr_2cols;
    my $my_filtered_gff         = $out_cds_intron_gff;
    ### ### [<now>] Finished extract_cds_introns...
    return [
            $my_cds_filtered_tbl,     $my_intron_filtered_tbl,
            $my_locusid_filtered_tbl, $my_filtered_gff,
            0
           ];

    #    }    #END ELSE IF NO SEQS  ARE INFRAME

}    #sub extract_cds_introns

## FUNCTION TO EXTRACT AND PROCESS SPLICE SITES AND START CODON
sub extract_sites ($my_input_nodots_gff, $locusid_2cols)
{

    ### [<now>] Running extract_sites at <file>[<line>]...

    my $my_command = "";
    ## SPLICE SITES
    #    print {*STDERR} "\nEXTRACT START AND SPLICE SITES from transcripts\n\n";

    #print {*STDERR} "$locus_id and $input_gff_fn\n";
    my @newsites = ();
    my $count    = 0;

    open(my $FH_LOC_sites, "<", "$locusid_2cols") or croak "Failed here";
    while (<$FH_LOC_sites>)
    {
        my ($id_genomic, $id_gene) = split;

        #  print {*STDERR} "$id_genomic,$id_gene\n";
        run("egrep -w '$id_gene\$' $my_input_nodots_gff > $tmp_dir/$id_gene.gff"
        );

        #  print {*STDERR} "$id_gene $input_gff_fn $tmp_dir/$id_gene.gff \n\n";
        ## POTENTIAL BUG SPLIT
        $my_command =
          "./bin/ssgff -dabeE $fastas_dir/$id_genomic $tmp_dir/$id_gene.gff > $tmp_dir/${id_gene}.all_sites";
        run($my_command);

        #~ foreach my $site (qw(Acceptor Donor Stop Start))
        #~ {

        #~ #    print {*STDERR} "egrep -A 1 $site $tmp_dir/${id_gene}.all_sites $sitesdir/${site}_sites.fa\n";
        #~ ## POTENTIAL BUG, split command below
        #~ $my_command = "grep -A 1 ";

        #~ run(" egrep -A 1 $site $tmp_dir/${id_gene}.all_sites | sed -e '/--/d' -e '/^\$/d' >> $sites_dir/${site}_sites.fa"
        #~ );
        #~ }
        $count++;
        print {*STDERR} "$count..";
    }    #while $FH_LOC_sites
    close $FH_LOC_sites;

    #    my $acceptors_fa  = "$sites_dir/Acceptor_sites.fa";
    #    my $acceptors_tbl = "$work_dir/Acceptor_sites.tbl";
    #
    #    my $donors_fa  = "$sites_dir/Donor_sites.fa";
    #    my $donors_tbl = "$work_dir/Donor_sites.tbl";
    #
    #    my $preATGx_fa  = "$sites_dir/Start_sites.fa";
    #    my $preATGx_tbl = "$work_dir/Start_sites.tbl";

    ## 2016.12.14a
    ## foreach my $site (qw(Acceptor Donor Start Stop ))
    foreach my $site (qw(Donor Acceptor  Start  ))
    {
        my $output_sites_fas = "$sites_dir/${site}_sites.fa";
        my $output_sites_tbl = "$work_dir/${site}_sites.tbl";
        $my_command =
          "grep -h -A 1 $site $tmp_dir/*.all_sites | grep -v '\-' > $output_sites_fas";
        run($my_command);
        ### [<now>] sites fasta_2_tbl at <file>[<line>]...
        fasta_2_tbl($output_sites_fas, $output_sites_tbl);
        if ($site eq 'Start')
        {
            my $ATGx_tbl = "$sites_dir" . "Start_sites_complete.tbl";
            $my_command =
              ##  "gawk '{printf \$1\" \";for (i=1;i<=60-length(\$2);i++) printf \"n\"; print \$2}' $output_sites_tbl  > $ATGx_tbl";
              "gawk '{printf \$1\" \";for (i=1;i<=60-length(\$2);i++) printf \"n\"; print \$2}' $output_sites_tbl  > $ATGx_tbl";
            run($my_command);
            $my_command =
              "cp $output_sites_tbl ATG_safe.tbl; mv $ATGx_tbl $output_sites_tbl ";
            run($my_command);
        }

        my ($canonical_tbl_fn, $noncanonical_num) =
          noncanonical($site, %canonical_const);

        #~ say "BAD GUYS";
        #~ print Dump @bad_guys;
        push @newsites, $canonical_tbl_fn;
        push @newsites, $noncanonical_num;

    }

    #~ push @newsites, $canonical_tbl_fn;
    #~ push @newsites, $noncanonical_num;

    #~ ##ADD N TO START SITES##
    #~ ## POTENTIAL BUG
    #~ my $ATGx_tbl = "$sites_dir" . "Start_sites_complete.tbl";
    #~ my $my_command =
    #~ "gawk '{printf \$1\" \";for (i=1;i<=60-length(\$2);i++) printf \"n\"; print \$2}' $preATGx_tbl > $ATGx_tbl";
    #~ run($my_command);

    #`gawk '{printf \$1" ";for (i=1;i<=60-length(\$2);i++) printf "n"; print \$2}' $prestarttbl > $sites_dir/Start_sites_complete.tbl`;

    #~ #

    #~ print {*STDERR} "\n\nEliminate non-canonical donors/acceptors/starts:\n";

    #~ #      ##EXTRACT NON CANONICAL DONORS
    #~ #my $noncanonical     = "";
    #~ my $generic_noncanonical       = "";
    #~ my $total_generic_noncanonical = "";
    #~ my $total_generic_canonical    = "";

    #~ my $new_donor_tbl = "";
    #~ $my_command =
    #~ "gawk '{print \$2}' $donors_tbl  | egrep -v '^[NATCGn]{31}GT' ";
    #~ my $noncanonical = capture($my_command);

    #~ my $donor_noncanonical_temp = $work_dir . $species . "_non_canonical_donor";
    #~ open(my $FH_FOUT, ">", "$donor_noncanonical_temp") or croak "Failed here";
    #~ print {$FH_FOUT} "$noncanonical";
    #~ close $FH_FOUT;

    #~ $total_generic_noncanonical =
    #~ num_of_lines_in_file($donor_noncanonical_temp);

    #~ print {*STDERR}
    #~ "\nThere are $total_generic_noncanonical non-canonical donors within the training set:\n";

    #~ #
    #~ if ($total_generic_noncanonical)
    #~ {    #if there are non canonical donors

    #~ my @generic_noncanonical = ();
    #~ open $FH_LOCID, "-|",
    #~ "egrep -wf $donor_noncanonical_temp $donors_tbl | gawk '{print \$1}' - | sort | uniq";
    #~ while (<$FH_LOCID>)
    #~ {

    #~ push(@generic_noncanonical, "$_");

    #~ }
    #~ close $FH_LOCID;

    #~ foreach my $line (@generic_noncanonical)
    #~ {

    #~ #   my (@noncanonical_array)= split (/\.\d+:/, $line);
    #~ my (@noncanonical_array) = split(/:/, $line);
    #~ my $first = $noncanonical_array[0] . ":";
    #~ $generic_noncanonical .= "$first\n";

    #~ }

    #~ #  unlink $donor_noncanonical_temp;

    #~ my $noncanonical_name_tmp =
    #~ $work_dir . $species . "_non_canonical_donor_seq_name";
    #~ open(my $FH_FOUT, ">", "$noncanonical_name_tmp")
    #~ or croak "Failed here";
    #~ print {$FH_FOUT} "$generic_noncanonical";
    #~ close $FH_FOUT;

    #~ open $FH_LOCID, "-|", "egrep -vf $noncanonical_name_tmp $donors_tbl";
    #~ while (<$FH_LOCID>)
    #~ {
    #~ $new_donor_tbl .= $_;
    #~ }
    #~ close $FH_LOCID;

    #~ my $tmp_canonical_donor_tbl_fn =
    #~ $work_dir . $species . ".canonical.donor.tbl";
    #~ open(my $FH_FOUT, ">", "$tmp_canonical_donor_tbl_fn")
    #~ or croak "Failed here";
    #~ print {$FH_FOUT} "$new_donor_tbl";
    #~ close $FH_FOUT;

    #~ # unlink $noncanonical_name_tmp;

    #~ $total_generic_canonical =
    #~ num_of_lines_in_file($tmp_canonical_donor_tbl_fn);

    #~ print {*STDERR}
    #~ "\nThere are $total_generic_canonical canonical donors within the training set:\n";

    #~ push(@newsites, "$tmp_canonical_donor_tbl_fn");
    #~ push(@newsites, "$total_generic_noncanonical");
    #~ }
    #~ else
    #~ {    #if there are no non-canonical
    #~ my $total_generic_canonical = num_of_lines_in_file($donors_tbl);

    #~ print {*STDERR}
    #~ "There are $total_generic_canonical canonical donors within the training set:\n";
    #~ push(@newsites, "$donors_tbl");
    #~ push(@newsites, "");

    #~ }    #if there are no non-canonical

    #~ #      #

    #~ #      ##
    #~ #      ##EXTRACT NON CANONICAL ACCEPTORS
    #~ $noncanonical            = "";
    #~ $generic_noncanonical    = "";
    #~ $total_generic_canonical = "";
    #~ my $acceptor_new_tbl = "";

    #~ my $foobar_tmp = "";

    #~ #$my_command =
    #~ #  "gawk '{print \$2}' $acceptortbl | egrep -v '^[NATCG]{28}AG'";

    #~ # BUG this blows if there are no such sites...
    #~ #$foobar_tmp = capture($my_command);

    #~ open $FH_LOCID, "-|",
    #~ "gawk '{print \$2}' $acceptors_tbl | egrep -v '^[NATCG]{28}AG'";
    #~ while (<$FH_LOCID>)
    #~ {
    #~ $noncanonical .= $_;
    #~ }
    #~ if (length($noncanonical) > 0)
    #~ {
    #~ close $FH_LOCID;
    #~ }

    #~ my $acceptor_noncanonical_temp =
    #~ $work_dir . $species . "_non_canonical_acceptor";
    #~ open(my $FH_FOUT, ">", "$acceptor_noncanonical_temp") or croak "Failed here";
    #~ print {$FH_FOUT} "$noncanonical";
    #~ close $FH_FOUT;

    #~ $total_generic_noncanonical =
    #~ num_of_lines_in_file($acceptor_noncanonical_temp);

    #~ print {*STDERR}
    #~ "\nThere are $total_generic_noncanonical non-canonical acceptors within the training set:\n";
    #~ #
    #~ if ($total_generic_noncanonical)
    #~ {    #if there are non-canonical acceptors

    #~ my @generic_noncanonical = ();
    #~ open $FH_LOCID, "-|",
    #~ "egrep -f $acceptor_noncanonical_temp $acceptors_tbl | gawk '{print \$1}' - | sort | uniq ";
    #~ while (<$FH_LOCID>)
    #~ {
    #~ push(@generic_noncanonical, "$_");
    #~ }

    #~ close $FH_LOCID;

    #~ foreach my $line (@generic_noncanonical)
    #~ {
    #~ my (@noncanonical_array) = split(/:/, $line);
    #~ my $first = $noncanonical_array[0] . ":";
    #~ $generic_noncanonical .= "$first\n";

    #~ }

    #~ unlink $acceptor_noncanonical_temp;

    #~ my $noncanonical_name_tmp =
    #~ $work_dir . $species . "_non_canonical_acceptor_seq_name";

    #~ open(my $FH_FOUT, ">", "$noncanonical_name_tmp");
    #~ print {$FH_FOUT} "$generic_noncanonical";
    #~ close $FH_FOUT;

    #~ open $FH_LOCID, "-|", "egrep -vf $noncanonical_name_tmp $acceptors_tbl";
    #~ while (<$FH_LOCID>)
    #~ {
    #~ $acceptor_new_tbl .= $_;
    #~ }
    #~ close $FH_LOCID;

    #~ #unlink $noncanonical_name_tmp;

    #~ my $acceptor_canonical_temp =
    #~ $work_dir . $species . ".canonical.acceptor.tbl";
    #~ open(my $FH_FOUT, ">", "$acceptor_canonical_temp")
    #~ or croak "Failed here";
    #~ print {$FH_FOUT} "$acceptor_new_tbl";
    #~ close $FH_FOUT;

    #~ #unlink $noncanonical_name_tmp;

    #~ my $total_generic_canonical =
    #~ num_of_lines_in_file($acceptor_canonical_temp);

    #~ print {*STDERR}
    #~ "\nThere are $total_generic_canonical canonical acceptors within the training set:\n";

    #~ push(@newsites, "$acceptor_canonical_temp");
    #~ push(@newsites, "$total_generic_noncanonical");

    #~ }
    #~ else
    #~ {    #if there are only canonical use initial file list
    #~ my $total_generic_canonical = num_of_lines_in_file($acceptors_tbl);

    #~ print {*STDERR}
    #~ "There are $total_generic_canonical canonical acceptors within the training set:\n";
    #~ push(@newsites, "$acceptors_tbl");
    #~ push(@newsites, "0");
    #~ }    #if there are only canonical use initial file list

    #~ #      #

    #~ #      #
    #~ #      ##EXTRACT NON CANONICAL STARTS

    #~ $noncanonical            = "";
    #~ $generic_noncanonical    = "";
    #~ $total_generic_canonical = "";
    #~ my $new_start_tbl = "";

    #~ open $FH_LOCID, "-|",
    #~ "gawk '{print \$2}' $ATGx_tbl | egrep -v '^[NATCG]{30}ATG' ";
    #~ while (<$FH_LOCID>)
    #~ {
    #~ $noncanonical .= $_;
    #~ }
    #~ if (length($noncanonical) > 0)
    #~ {
    #~ close $FH_LOCID;
    #~ }

    #~ #close $FH_LOCID;

    #~ my $temp_startnoncanonical = $work_dir . $species . "_non_canonical_start";
    #~ open(my $FH_FOUT, ">", "$temp_startnoncanonical") or croak "Failed here";
    #~ print {$FH_FOUT} "$noncanonical";
    #~ close $FH_FOUT;

    #~ $total_generic_noncanonical = num_of_lines_in_file($temp_startnoncanonical);

    #~ print {*STDERR}
    #~ "\nThere are $total_generic_noncanonical non-canonical starts within the training set:\n";
    #~ #

    #~ if ($total_generic_noncanonical)
    #~ {    #if there are non-canonical starts

    #~ my @generic_noncanonical = ();
    #~ open $FH_LOCID, "-|",
    #~ "egrep -wf $temp_startnoncanonical $ATGx_tbl | gawk '{print \$1}' - | sort | uniq ";
    #~ while (<$FH_LOCID>)
    #~ {
    #~ push(@generic_noncanonical, "$_");
    #~ }
    #~ close $FH_LOCID;

    #~ foreach my $line (@generic_noncanonical)
    #~ {
    #~ my (@noncanonical_array) = split(/:/, $line);
    #~ my $first = $noncanonical_array[0] . ":";
    #~ $generic_noncanonical .= "$first\n";

    #~ }

    #~ unlink $temp_startnoncanonical;

    #~ my $noncanonical_name_tmp =
    #~ $work_dir . $species . "_non_canonical_start_seq_name";
    #~ open(my $FH_FOUT, ">", "$noncanonical_name_tmp");
    #~ print {$FH_FOUT} "$generic_noncanonical";
    #~ close $FH_FOUT;

    #~ open(my $FH_LOCID, "-|", "egrep -vf $noncanonical_name_tmp $ATGx_tbl ");
    #~ while (<$FH_LOCID>)
    #~ {
    #~ $new_start_tbl .= $_;
    #~ }
    #~ close $FH_LOCID;

    #~ # unlink $noncanonical_name_tmp;

    #~ my $tmp_canonical_ATGx_tbl_fn =
    #~ $work_dir . $species . ".canonical.start.tbl";
    #~ open(my $FH_FOUT, ">", "$tmp_canonical_ATGx_tbl_fn")
    #~ or croak "Failed here";
    #~ print {$FH_FOUT} "$new_start_tbl";
    #~ close $FH_FOUT;

    #~ #unlink $noncanonical_name_tmp;
    #~ my $total_generic_canonical =
    #~ num_of_lines_in_file($tmp_canonical_ATGx_tbl_fn);

    #~ print {*STDERR}
    #~ "\nThere are $total_generic_canonical canonical starts within the training set:\n";

    #~ push(@newsites, "$tmp_canonical_ATGx_tbl_fn");
    #~ push(@newsites, "$total_generic_noncanonical");

    #~ }
    #~ else
    #~ {
    #~ my $total_generic_canonical = num_of_lines_in_file($ATGx_tbl);

    #~ print {*STDERR}
    #~ "\nThere are $total_generic_canonical canonical starts within the training set:\n";
    #~ push(@newsites, "$ATGx_tbl");
    #~ push(@newsites, "0");

    #~ }    #if there are only canonical starts
    #~ #
    ## DEBUG
    #~ say "NEWSITES START";
    #~ print Dump @newsites;
    #~ say "NEWSITES END";
    ### [<now>] Finished extract_sites...
    return \@newsites;

}    #subectractprocesssites end

sub noncanonical ($site_type, %canonical_const)
{
    ### [<now>] Running noncanonical at <file>[<line>]...
    my $input_site_tbl   = "$work_dir/${site_type}_sites.tbl";
    my $out_site_tbl     = "$work_dir/${site_type}_canonical.DK.tbl";
    my $pre_bases_number = $canonical_const{$site_type}[0];
    my $pattern          = $canonical_const{$site_type}[1];
    my $pattern_len      = length($pattern);
    ## DEBUG
    #~ print Dump $pre_bases_number;
    #~ print Dump $pattern;
    #~ print Dump $pattern_len;
    my $canonical_counter     = 0;
    my $non_canonical_counter = 0;
    my @non_canonical_list    = ();

    open(my $FH_TBL_IN,  "<", "$input_site_tbl") or croak "Failed here";
    open(my $FH_TBL_OUT, ">", "$out_site_tbl")   or croak "Failed here";
    while (<$FH_TBL_IN>)
    {
        my ($site_id, $site_whole_seq) = split;
        my $test_seq = substr($site_whole_seq, $pre_bases_number, $pattern_len);
        ## DEBUG
        ## say $site_whole_seq;
        ## say $test_seq;
        if ($test_seq eq $pattern)
        {
            print {$FH_TBL_OUT} "$site_id\t$site_whole_seq\n";
            $canonical_counter += 1;
        }
        else
        {
            ## model10054m000225.10:10054.m000225
            push @non_canonical_list, $site_id;
            $non_canonical_counter += 1;
        }
    }

    close $FH_TBL_IN;
    close $FH_TBL_OUT;
    ### [<now>] Finished noncanonical...
    #return $canonical_counter, $non_canonical_counter; #, \@non_canonical_list;
    return $out_site_tbl, $non_canonical_counter;
}

## FUNCTION TO OBTAIN MARKOV MODELS CORRESPONDING TO THE CODING POTENTIAL
sub derive_coding_potential ($in_cds_tbl_fn, $in_intron_tbl_fn)
{
    ### [<now>] Running derive_coding_potential at <file>[<line>]...

    my $markov_mod_A = "";
    my $markov_mod_B = "";

    my $my_command =
      "gawk '{ l=length(\$2); L+=l;} END{ print L;}' $in_cds_tbl_fn";
    my $total_codingbases = capture($my_command);
    chomp $total_codingbases;

    $my_command =
      "gawk '{ l=length(\$2); L+=l;} END{ print L;}' $in_intron_tbl_fn ";
    my $total_noncodingbases = capture($my_command);
    chomp $total_noncodingbases;

    print {*STDERR}
      "There are $total_codingbases coding bases and $total_noncodingbases non-coding bases on this training set:\n";

    if (
        (
            $total_codingbases > $coding_bp_limit_A
         && $total_noncodingbases > $non_coding_bp_limit_A
        )
        || (   $total_codingbases > $coding_bp_limit_B
            && $total_noncodingbases > $non_coding_bp_limit_B)
        || (   $total_noncodingbases > $non_coding_bp_limit_C
            && $total_codingbases >
            ($multi_total_noncodingbases * $total_noncodingbases))
       )
    {
        $markov_mod_A = 5;
        $markov_mod_B = 4;
        print {*STDERR}
          "Deriving a markov model of order $markov_mod_A OPTION_1\n";

    }
    else
    {
        $markov_mod_A = 5;
        $markov_mod_B = 4;
        print {*STDERR}
          "Deriving a markov model of order $markov_mod_A  OPTION_2\n";
    }

    open(my $FH_INTRONS_tbl, "<", "$in_intron_tbl_fn") or croak "Failed here";
    my @intron_seqs = ();
    while (<$FH_INTRONS_tbl>)
    {
        my @columns_i = split;

        #print {*STDERR} "SECOND FIELD $i[2]";
        push @intron_seqs, $columns_i[1];
    }
    close $FH_INTRONS_tbl;

    open(my $FH_CDSes_tbl, "<", "$in_cds_tbl_fn") or croak "Failed here";
    my @coding_seqs = ();
    while (<$FH_CDSes_tbl>)
    {
        my @columns = split;
        push @coding_seqs, $columns[1];
    }
    close $FH_CDSes_tbl;

    print {*STDERR} "Intron model\n markov: ($markov_mod_B)";

    my $intron_initial =
      geneidCEGMA::SequenceModel->new('intron', 'FREQ', $markov_mod_B,
                                      \@intron_seqs, $intron_param_X, 0);

    my $intron_transition =
      geneidCEGMA::SequenceModel->new('intron', 'MM', $markov_mod_A,
                                      \@intron_seqs, $intron_param_X, 0);

    print {*STDERR} "Coding model\n";

    my $coding_initial =
      geneidCEGMA::SequenceModel->new('coding', 'FREQ', $markov_mod_A - 1,
                                      \@coding_seqs, $coding_param_X, 2);

    my $coding_transition =
      geneidCEGMA::SequenceModel->new('coding', 'MM', $markov_mod_A,
                                      \@coding_seqs, $coding_param_X, 2);

    my $initial_logs = geneidCEGMA::log_ratio($coding_initial, $intron_initial);

    my $transition_logs =
      geneidCEGMA::log_ratio($coding_transition, $intron_transition);

    geneidCEGMA::write_log($initial_logs, "$work_dir/coding.initial.5.logs");

    geneidCEGMA::write_log($transition_logs,
                           "$work_dir/coding.transition.5.logs");

    open(my $FH_PROFILE_1, "<", "$work_dir/coding.initial.5.logs")
      or croak "Failed here";
    my @profile_init = ();
    while (<$FH_PROFILE_1>)
    {
        last if m/^\s/;
        last if m/^[^ACGT]/;    #no actg
                                #last if m/^[^ACGTacgt]/;
        next if m/^#/;
        chomp;
        my @columns_profile = split;
        push @profile_init, \@columns_profile;
    }
    close $FH_PROFILE_1;

    open(my $FH_PROFILE_2, "<", "$work_dir/coding.transition.5.logs")
      or croak "Failed here";
    my @profile_trans = ();
    ## TODO code dupl 1
    while (<$FH_PROFILE_2>)
    {
        last if m/^\s/;
        last if m/^[^ACGT]/;    #no actg
                                #last if m/^[^ACGTacgt]/;
        next if m/^#/;
        chomp;
        my @columns_profile = split;
        push @profile_trans, \@columns_profile;
    }
    close $FH_PROFILE_2;
    ### [<now>] Finished derive_coding_potential...
    return [
            \@profile_init,        \@profile_trans, $total_codingbases,
            $total_noncodingbases, $markov_mod_A
           ];

}    #derive coding potential

## PROCESS SEQUENCES FUNCTION ( FLANKED GENE MODELS OBTAINED FOR OPTIMIZATION)
sub process_seqs_4opty ($my_input_nodots_gff, $opt_type, $contig_opt_flag)
{

    ### [<now>] Running process_seqs_4opty at <file>[<line>]...

    #my $pso_out_tbl = "";
    my $pso_gp_tbl  = "";
    my $gp_from_gff = "";

    #~ my $gp_fasta    = ""; #unused
    my $pso_gp_gff = "";
    my $my_command = "";

    #my $work_dir;

    open(
         my $FH_LOCID, "-|",
         "./bin/gff2gp.pypy $my_input_nodots_gff | sort -k 1 " );
         ## 2017.07.11
         ##"./bin/gff2gp.awk $my_input_nodots_gff | sort -k 1 "
        ##);

    while (<$FH_LOCID>)
    {

        $gp_from_gff .= $_;
    }
    close $FH_LOCID;

    my $pso_tmp_gp_from_gff =
      $work_dir . $species . $opt_type . $contig_opt_flag . ".gp";

    open(my $FH_FOUT, ">", "$pso_tmp_gp_from_gff");
    print {$FH_FOUT} "$gp_from_gff";
    close $FH_FOUT;

    # print {*STDERR} "BEFORE GETGENES: $fastas_dir, $pso_tmp_gp_from_gff, $work_dir/, $pso_out_tbl\n";
    print {*STDERR}
      "BEFORE GETGENES: $fastas_dir, $pso_tmp_gp_from_gff, $work_dir , $contig_opt_flag \n";

    my $gp_Xgetgenes_tmp_pre_tbl =
      get_genes($fastas_dir, $pso_tmp_gp_from_gff, $work_dir);
    print {*STDERR} "PRETBL AFTER GETGENES: $gp_Xgetgenes_tmp_pre_tbl \n";

    print {*STDERR}
      "\nGet sequences of 400-nt flanked sequences in tabular and gff formats\n";

    #~ my $seq4Optimization_temp_1_fn =  "$work_dir/processSequences4Optimization_temp1.txt";
    #~ $my_command =  "gawk 'BEGIN{b=\"x\"}{if (\$1!=b){printf \"\\n\%s \",\$1}printf \"\%s\",\$6;b=\$1}END{printf \"\\n\"}' $gp_Xgetgenes_tmp_pre_tbl > $seq4Optimization_temp_1_fn";
    #~ run($my_command);
    #~ $my_command =  "sort seq4Optimization_temp_1_fn | sed '/^\$/d' - | gawk '{print \$1, toupper(\$2)}' - |";
    #open( $FH_LOCID, $my_command );
    open(
        $FH_LOCID,
        "-|",
        "gawk 'BEGIN{b=\"x\"}{if (\$1!=b){printf \"\\n\%s \",\$1}printf \"\%s\",\$6;b=\$1}END{printf \"\\n\"}' $gp_Xgetgenes_tmp_pre_tbl | sort | sed '/^\$/d' - | gawk '{print \$1, toupper(\$2)}' - "
    );

    while (<$FH_LOCID>)
    {
        $pso_gp_tbl .= $_;
    }
    close $FH_LOCID;

    my $temp_gp_tbl =
      $work_dir . $species . $opt_type . $contig_opt_flag . ".gp.tbl";
    open($FH_FOUT, ">", "$temp_gp_tbl") or croak "Failed here";
    print {$FH_FOUT} "$pso_gp_tbl";
    close $FH_FOUT;

    open(
        $FH_LOCID,
        "-|",
        "gawk 'BEGIN{OFS=\"\\t\";pos=1;b=\"x\"}{if (\$1!=b){pos=1}; print \$1,\"annotations\",\$3,pos,pos+\$5-1,\"\.\",\"+\",\"\.\",\$1\$2; pos+=\$5;b=\$1 }' $gp_Xgetgenes_tmp_pre_tbl | egrep -v '(Intron|Utr)' - "
    );
    while (<$FH_LOCID>)
    {
        $pso_gp_gff .= $_;
    }
    close $FH_LOCID;

    my $temp_gp_gff =
      $work_dir . $species . $opt_type . $contig_opt_flag . ".gp.gff";
    open($FH_FOUT, ">", "$temp_gp_gff") or croak "Failed here";
    print {$FH_FOUT} "$pso_gp_gff";
    close $FH_FOUT;

    print {*STDERR} "DONE\n";

    print {*STDERR}
      "\nGet sequences of 400-nt flanked sequences in multi-fasta format\n";

    my $temp_gp_fa =
      $work_dir . $species . $opt_type . $contig_opt_flag . ".gp.fa";

    $temp_gp_fa = tbl_2_fasta($temp_gp_tbl, $temp_gp_fa);

    ## BUG keep intermediates
    #~ unlink $gp_Xgetgenes_tmp_pre_tbl;

    print {*STDERR} "\nSet up files for optimization\n\n";
    $my_command = "gawk '{print \$1,length(\$2)}' $temp_gp_tbl | sort -k1,1 ";
    my $gp_seqs_lengths = capture($my_command);

    #  ` gawk '{print \$1,length(\$2)}' $temp_gp_tbl | sort -k1,1 `;    ##XX

    my $opt_type_cds_contigs_total_bp =
      $work_dir . $species . $opt_type . $contig_opt_flag . ".gp_cds_length";
    open($FH_FOUT, ">", "$opt_type_cds_contigs_total_bp")
      or croak "Failed here";
    print {$FH_FOUT} "$gp_seqs_lengths";
    close $FH_FOUT;

    my $cds_gp = "";
    open(
        $FH_LOCID,
        "-|",
        "./bin/gff2cds.awk source=\"annotations\" $temp_gp_gff | sort -k1,1 | join $opt_type_cds_contigs_total_bp - "
    );
    while (<$FH_LOCID>)
    {
        $cds_gp .= $_;
    }
    close $FH_LOCID;

    my $temp_gp_cds_fn =
      $work_dir . $species . $opt_type . $contig_opt_flag . ".cds_gp";
    open($FH_FOUT, ">", "$temp_gp_cds_fn") or croak "Failed here";
    print {$FH_FOUT} "$cds_gp";
    close $FH_FOUT;

    ## 2016.12.14a unused ???
    ##my $pso_gp_eval_gff = "";

    my $eval_gp_out_gff =
      $work_dir . $species . $opt_type . $contig_opt_flag . ".gp_eval_gff";
    $my_command =
      "./bin/1800p_script.awk $temp_gp_cds_fn $temp_gp_gff > $eval_gp_out_gff";
    run($my_command);

    #open( $FH_LOCID, "-|",
    # "gawk 'BEGIN{while (getline<ARGV[1]>0){len[\$1]=\$2;};ARGV[1]=\"\";OFS=\"\\t\";}{if (NR==1) {ant=\$1;print \$1,\$2,\"Sequence\",1,len[ant],\"\.\",\"\.\",\"\.\",\"\.\"};if (\$1!=ant) {print \"\#\$\";ant=\$1;print \$1,\$2,\"Sequence\",1,len[ant],\"\.\",\"\.\",\"\.\",\"\.\"}; print }' $temp_gp_cds_fn $temp_gp_gff "
    #    );
    #while (<$FH_LOCID>) {
    #    $pso_gp_eval_gff .= $_;

    #    }
    #    close $FH_LOCID;

    #    my $eval_gp_out_gff = $work_dir . $species . $opt_type . ".gp_eval_gff";
    #    open( $FH_FOUT, ">", "$eval_gp_out_gff" ) or croak "Failed here";
    #    print {$FH_FOUT} "$pso_gp_eval_gff";
    #    close $FH_FOUT;

    ##if ($contig_opt_flag)
    ##{

    #my $temp_gp_fa = $species.$opt_type.".gp.fa";

    #$temp_gp_fa = FastaToTbl($temp_gp_tbl,$temp_gp_fa);

    my @gp_tabular = split(/\n/, $pso_gp_tbl);
    my $seq = "";
    foreach my $line (@gp_tabular)
    {
        #        chomp $line;
        #        my @columns_f = split " ", $line;

        my @columns_f = split /\s+/, $line;
        $seq .= $columns_f[1];
    }
    my $my_seq_len    = length($seq);
    my $folded_seq_gp = fold4fasta($seq);
    my $temp_fastagpcontig =
      $work_dir . $species . $opt_type . $contig_opt_flag . ".combined.gp.fa";
    open($FH_FOUT, ">", "$temp_fastagpcontig") or croak "Failed here";
    print {$FH_FOUT} ">$species\n$folded_seq_gp\n";
    close $FH_FOUT;

    my $temp_tabulargpcontig =
      $work_dir . $species . $opt_type . $contig_opt_flag . ".combined.gp.tbl";
    open($FH_FOUT, ">", "$temp_tabulargpcontig") or croak "Failed here";
    print {$FH_FOUT} "$species\t$seq\n";
    close $FH_FOUT;
    $my_command = "gawk '{print \$1,length(\$2)}' $temp_tabulargpcontig";
    my $gp_contigs_seqs_lengths = capture($my_command);

    #` gawk '{print \$1,length(\$2)}' $temp_tabulargpcontig `;

    my $temp_seqlencontig =
        $work_dir
      . $species
      . $opt_type
      . $contig_opt_flag
      . ".gp_cds_contig_length";
    open($FH_FOUT, ">", "$temp_seqlencontig") or croak "Failed here";
    print {$FH_FOUT} "$gp_contigs_seqs_lengths";
    close $FH_FOUT;

    my $gp_contig_tmp = "";
    open(
        $FH_LOCID,
        "-|",
        "./bin/multiple_annot2one.awk species=$species leng=$my_seq_len $temp_gp_cds_fn "
    );
    while (<$FH_LOCID>)
    {

        $gp_contig_tmp .= $_;
    }
    close $FH_LOCID;

    my $temp_gff2gpcontig =
      $work_dir . $species . $opt_type . $contig_opt_flag . ".contig.gp.cds";

    open($FH_FOUT, ">", "$temp_gff2gpcontig") or croak "Failed here";
    print {$FH_FOUT} "$gp_contig_tmp";
    close $FH_FOUT;

    my $cds2gffcontig = "";
    open(
        $FH_LOCID,
        "-|",
        "./bin/cds2gff.awk $temp_gff2gpcontig | gawk 'BEGIN{OFS=\"\\t\";}{if (NR==1){print \"$species\",\"annotations\",\"Sequence\",\"1\",$my_seq_len,\".\",\".\",\".\",\".\";print}else {print}}' - "
    );
    while (<$FH_LOCID>)
    {
        $cds2gffcontig .= $_;
    }
    close $FH_LOCID;

    my $temp_gp_cdsgff_contig_eval =
        $work_dir
      . $species
      . $opt_type
      . $contig_opt_flag
      . ".cds_gp_contig.eval.gff";
    open($FH_FOUT, ">", "$temp_gp_cdsgff_contig_eval")
      or croak "Failed here";
    print {$FH_FOUT} "$cds2gffcontig";
    close $FH_FOUT;
    ## TODO Inspect files / values returned. 2016.12.15
    ### Finished process_seqs_4opty...
    return [
            $temp_gp_cdsgff_contig_eval, $temp_fastagpcontig,
            $temp_tabulargpcontig,       $temp_seqlencontig
           ];

    ##}
    #    elsif (!$contig_opt_flag)
    #    {
    #        print {*STDERR} "L1803, NOT CONTIG OPT\n";
    #        return [
    #            $eval_gp_out_gff, $temp_gp_fa,
    #            $temp_gp_tbl,     $opt_type_cds_contigs_total_bp
    #        ];
    #    }

}    #processSequences optimization

## GETGENES FUNCTION: EXTRACT FLANKED SEQUENCES FROM GENE MODELS FOR LATER OPTIMIZATION
sub get_genes ($my_fastas_dir, $my_pso_tmp_gp_from_gff, $work_dir)
{
    ### [<now>] Running get_genes at <file>[<line>]...
    my $only_non_red = 0;
    ## unused
    ## my $prev_gene   = "x";
    ## my $prev_chrom   = "x";
    ## my $trail      = "";

    #my %genenames; unused var
    ## DEBUG: testing chrom fasta_2_tbl
    my $gp_out_tbl = "$work_dir/gp_exon_utr_intron.tbl";

    #~ my $chrom_tmp_tbl = "$tmp_dir/chrom_tmp.tbl";
    #~ ### [<now>] chrom fasta_2_tbl at <file>[<line>]...
    #~ $chrom_tmp_tbl =
    #~ fasta_2_tbl($my_fastas_dir . $chrom, $chrom_tmp_tbl);
    ## DEBUG END...
    #~ chomp($my_fastas_dir);
    #~ chomp($genes_fn_X);
    #~ chomp($work_dir);
    #~ chomp($gp_out_tbl);

    open(my $FH_REFGENE, "<", "$my_pso_tmp_gp_from_gff")
      or croak "Failed here";
    ### [<now>] Reading at <file>[<line>]...
    ###  $my_pso_tmp_gp_from_gff

    open(my $FH_OUT_tblgb, ">", "$gp_out_tbl") or croak "Failed here";
    while (<$FH_REFGENE>)
    {

        #my $line = <$FH_REFGENE>;
        #chomp $line;
        my @tmp_cols = split;    #, $_;
        ## DEBUG
        #~ print Dump @tmp_cols;
        my (
            $name,     $chrom,          $strand,
            $tx_start, $tx_end,         $cds_start,
            $cds_end,  $exon_count_int, $tmp_exons_starts,
            $tmp_exons_ends
           )
          = @tmp_cols;

        #= split " ", $line;

        my @exon       = split ",", $tmp_exons_starts;
        my @exons_ends = split ",", $tmp_exons_ends;

        ##( $tmp_exon =~ m/(\d+)/g );
        push @exon, @exons_ends;

        ## DEBUG
        #~ print Dump @exon;

        ##= split " ", $line;
        ##while (<$FH_REFGENE>) {
        #~ ##split
        #~ m/([\w\-\.:]+)\s+([\w\.\-:]+)\s+([\+\-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([^\n]+)/;

        #~ my $name          = $1;
        #~ my $chrom         = $2;
        #~ my $strand        = $3;
        #~ my $tx_start      = $4;
        #~ my $tx_end        = $5;
        #~ my $cds_start     = $6;
        #~ my $cds_end       = $7;
        #~ my $exon_count_int = $8;
        #~ my @exon          = ( $9 =~ m/(\d+)/g );

        my $cds_len    = $cds_end - $cds_start;
        my $tx_len     = $tx_end - $tx_start;
        my $cds_offset = $cds_start - $tx_start;

        my $redundant = 0;
        ## 2016.12.13c unused???
        ## my $i                = 0;    # exon counter
        my $j    = 0;
        my $call = "";
        ## 2016.12.14a unused here??
        ## my $sub_seq          = "";
        my $genomic_tmp_seqX = "";

        #my @tabular = ();

        if (!$only_non_red || ($only_non_red && !$redundant))
        {
            open(my $FH_FLEN, "<", "$my_fastas_dir${chrom}_len")
              or croak "Failed here";

            ## DEBUG: fix the loop/while etc
            my $line = <$FH_FLEN>;
            chomp $line;
            my @le = split " ", $line;
            close $FH_FLEN;

            #added code
            my $chrom_tmp_tbl = "$tmp_dir/chrom_tmp.tbl";
            ### [<now>] chrom fasta_2_tbl at <file>[<line>]...
            $chrom_tmp_tbl =
              fasta_2_tbl($my_fastas_dir . $chrom, $chrom_tmp_tbl);

            #print {*STDERR} "FATOTBL: $my_fastas_dir"."$chro\n";
            open(my $FH_IN, "<", "$chrom_tmp_tbl") or croak "Failed here";

            #            my @tabular = ();
            my $sub_seq = "";
            while (<$FH_IN>)
            {
                my @columns_f = split;
                $sub_seq .= $columns_f[1];

                #chomp;
                #print {*STDERR} "$_";
                #                push @tabular, "$_";
            }
            close $FH_IN;

            #  print {*STDERR} "\nGP: @tabular\n";

            #            my $sub_seq = "";

            #my $sublen = 0;
            #            foreach my $line (@tabular)
            #            {
            #                chomp $line;
            #                my @columns_f = split " ", $line;
            #
            #                #print {*STDERR} "$f[0]\n";
            #                $sub_seq .= $columns_f[1];
            #
            #            }
            ## DEBUG added

            if ($le[1] < $tx_end)
            {

                my $new_len = $le[1];
                $genomic_tmp_seqX =
                  substr($sub_seq, $tx_start, ($new_len - $tx_start));

            }
            elsif ($le[1] >= $tx_end)
            {
                $genomic_tmp_seqX = substr($sub_seq, $tx_start, $tx_len);

            }

            # my $genomic_tmp_seqX = `$call`;

            my $genomic_len = length($genomic_tmp_seqX);
            my $cds_seq     = "";

            if ($genomic_len == 0)
            {
                print {*STDERR} "getgenes: gene of 0 length ($name), $call\n";
                next;
            }

            #  if ($genomic_len != $new_len) {
            #    print {*STDERR} "getgenes: length mismatch ($name)\n";
            #    next;
            #  }

            ##for (my $ii = 0 ; $ii < $exon_count_int ; $ii++)
            ## XXX fixing C style loop XXX
            ##my $ii = 0;
            for my $ii (0 .. $exon_count_int - 1)
            {
                my $utr_A           = 0;
                my $utr_B           = 0;
                my $utr_S           = 0;
                my $utr_L           = 0;
                my $exSt            = $exon[$ii] - $cds_start;
                my $exon_length_int = $exon[$ii + $exon_count_int] - $exon[$ii];
                my $ex_type         = "Internal";

                if ($exSt + $exon_length_int > 0 && $exSt < $cds_len)
                {    # cds

                    if ($exSt <= 0 || $ii == 0)
                    {
                        if ($strand eq '+')
                        {
                            $ex_type = "First";
                        }
                        else
                        {
                            $ex_type = "Terminal";
                        }
                    }

                    if (   $exSt + $exon_length_int >= $cds_len
                        || $ii == $exon_count_int - 1)
                    {
                        if ($strand eq '+')
                        {
                            $ex_type = "Terminal";
                        }
                        else
                        {
                            $ex_type = "First";
                        }
                    }

                    if ($exSt <= 0 && $exSt + $exon_length_int >= $cds_len)
                    {
                        $ex_type = "Single";
                    }

                    if ($exSt < 0)
                    {
                        $utr_B           = 1;
                        $utr_S           = $exSt;
                        $utr_L           = abs($exSt);
                        $exon_length_int = $exon_length_int - abs($exSt);
                        $exSt            = 0;
                    }

                    if ($exSt + $exon_length_int > $cds_len)
                    {
                        $utr_A = 1;
                        $utr_S = $cds_len;
                        $utr_L = $exon_length_int - ($cds_len - $exSt);
                        $exon_length_int = $cds_len - $exSt;
                    }

                    my $iex;
                    my $seq = substr($genomic_tmp_seqX, $exSt + $cds_offset,
                                     $exon_length_int);

                    $seq = lc($seq);

                    if ($strand eq '+')
                    {    # forward

                        if ($utr_B)
                        {
                            my $iutr = $ii + 1;
                            my $my_utrsX =
                              substr($genomic_tmp_seqX, $utr_S + $cds_offset,
                                     $utr_L);

                            $my_utrsX = lc($my_utrsX);
                            $cds_seq  = $cds_seq
                              . "$name\t$chrom\tUtr\t$iutr\t$utr_L\t$my_utrsX\n";
                        }

                        $iex     = $ii + 1;
                        $cds_seq = $cds_seq
                          . "$name\t$chrom\t$ex_type\t$iex\t$exon_length_int\t$seq\t$exon[$ii]\t$exon[$ii+$exon_count_int]\n";

                        if ($utr_A)
                        {
                            my $iutr = $ii + 1;
                            my $my_utrsX =
                              substr($genomic_tmp_seqX, $utr_S + $cds_offset,
                                     $utr_L);

                            $my_utrsX = lc($my_utrsX);
                            $cds_seq  = $cds_seq
                              . "$name\t$chrom\tUtr\t$iutr\t$utr_L\t$my_utrsX\n";
                        }

                    }
                    else
                    {    # reverse

                        if ($utr_B)
                        {
                            my $iutr = $exon_count_int - $ii;
                            my $my_utrsX =
                              substr($genomic_tmp_seqX, $utr_S + $cds_offset,
                                     $utr_L);

                            $my_utrsX = lc($my_utrsX);
                            $my_utrsX =~ tr/acgt/tgca/;
                            $my_utrsX = reverse($my_utrsX);
                            $cds_seq =
                              "$name\t$chrom\tUtr\t$iutr\t$utr_L\t$my_utrsX\n"
                              . $cds_seq;
                        }

                        $iex = $exon_count_int - $ii;
                        $seq =~ tr/acgt/tgca/;
                        $seq = reverse($seq);
                        $cds_seq =
                          "$name\t$chrom\t$ex_type\t$iex\t$exon_length_int\t$seq\t$exon[$ii+$exon_count_int]\t$exon[$ii]\n"
                          . $cds_seq;

                        if ($utr_A)
                        {
                            my $iutr = $exon_count_int - $ii;
                            my $my_utrsX =
                              substr($genomic_tmp_seqX, $utr_S + $cds_offset,
                                     $utr_L);

                            $my_utrsX = lc($my_utrsX);
                            $my_utrsX =~ tr/acgt/tgca/;
                            $my_utrsX = reverse($my_utrsX);
                            $cds_seq =
                              "$name\t$chrom\tUtr\t$iutr\t$utr_L\t$my_utrsX\n"
                              . $cds_seq;
                        }

                    }

                    if (
                        $ex_type ne "Single"
                        && (   ($ex_type ne "Terminal" && $strand eq '+')
                            || ($ex_type ne "First" && $strand eq '-'))
                       )
                    {

                        my $inSt = $exon[$ii + $exon_count_int] - $cds_start;
                        my $inLe =
                          $exon[$ii + 1] - $exon[$ii + $exon_count_int];

                        if ($inSt + $inLe > 0 && $inSt < $cds_len)
                        {

                            if ($inSt < 0)
                            {
                                print "getgenes: intron out of range! (1)\n";
                                exit(1);
                            }

                            if ($inSt + $inLe > $cds_len)
                            {
                                print "getgenes: intron out of range! (2)\n";
                                exit(1);
                            }

                            $seq =
                              substr($genomic_tmp_seqX, $inSt + $cds_offset,
                                     $inLe);
                            $seq = "\L$seq";

                            my $iIn;

                            if ($strand eq '+')
                            {    # forward
                                $iIn     = $j + 1;
                                $cds_seq = $cds_seq
                                  . "$name\t$chrom\tIntron\t$iIn\t$inLe\t$seq\n";
                            }
                            else
                            {
                                $iIn = $exon_count_int - $j - 1;
                                $seq =~ tr/acgt/tgca/;
                                $seq = reverse($seq);
                                $cds_seq =
                                  "$name\t$chrom\tIntron\t$iIn\t$inLe\t$seq\n"
                                  . $cds_seq;
                            }

                        }
                        else
                        {
                            print "getgenes.pl: intron out of range! (3)\n";
                            if ($inSt + $inLe <= 0)
                            {
                                print "getgenes.pl: intron in 5' UTR\n";
                            }
                            else
                            {
                                print "getgenes.pl: intron in 3' UTR\n";
                            }

                            exit(1);
                        }

                        $j = $j + 1;
                    }

                }
                else
                {    # UTRs

                    $exSt = $exon[$ii] - $tx_start;
                    say $ii, " RRR ", $exon_count_int, " MMM ", $exon[$ii];
                    $exon_length_int =
                      $exon[$ii + $exon_count_int] - $exon[$ii];

                    my $my_utrsX =
                      substr($genomic_tmp_seqX, $exSt, $exon_length_int);

                    if ($strand eq '+')
                    {    # forward
                        my $iutr = $ii + 1;

                        $my_utrsX = lc($my_utrsX);
                        $cds_seq =
                          $cds_seq
                          . "$name\t$chrom\tUtr\t$iutr\t$exon_length_int\t$my_utrsX\n";
                    }
                    else
                    {    # reverse
                        my $iutr = $exon_count_int - $ii;

                        $my_utrsX = lc($my_utrsX);
                        $my_utrsX =~ tr/acgt/tgca/;
                        $my_utrsX = reverse($my_utrsX);
                        $cds_seq =
                          "$name\t$chrom\tUtr\t$iutr\t$exon_length_int\t$my_utrsX\n"
                          . $cds_seq;
                    }

                }
            }

            print {$FH_OUT_tblgb} $cds_seq;

        }
        elsif ($only_non_red)
        {
            print {*STDERR} "$name\n";
        }

    }

    close $FH_OUT_tblgb;
    close $FH_REFGENE;
    ### [<now>] Finished get_genes...
    #$gp_out_tbl = $gp_out_tbl;
    return $gp_out_tbl;

}    #getgenes                         at <file>[<line>]

## GET BACKGROUND SEQUENCES (ex. 62 Kmer) used to compute log likelihoods

sub bit_score_graph ($info_output, $info_thresh, $offset)
{
    ### [<now>] Running bit_score_graph at <file>[<line>]...
    print {*STDERR}
      "bitscoregraph input:  $info_output, $info_thresh, $offset\n";
    my @info = ();    ##($offset - 1, $offset + 1);
    Readonly::Scalar my $equal_sign_multi => 30;
    ## 2016.12.14a
    ## my @fields;

    open(my $FH_INFO, "<", "$info_output") or croak "Failed here";
    while (<$FH_INFO>)
    {    ## 2016.12.14a: these information files do not haveanythyng but numbers
        ## not needed
        #~ next if m/^#/;
        #~ last if m/^\s/;
        #~ last if m/^[^\d]/;
        #~ chomp;

        #print {*STDERR} "QQQ prefields: $_";

        my ($matrix_position, $info_score) = split;
        ## @fields = split;
        #~ printf {*STDERR} "%2s %2.2f %s",
        #~ ($fields[0], $fields[1], "=" x int($fields[1] * 30));
        #~ printf {*STDERR} "%2s %2.2f %s",
        #~ ($matrix_position, $info_score, "=" x int($info_score * 30));

        printf {*STDERR} "%2s %2.2f %s",
          (
            $matrix_position, $info_score,
            "=" x int($info_score * $equal_sign_multi)
          );

        if ($info_score > $info_thresh)
        {
            push(@info, $matrix_position);
        }
        print {*STDERR} "\n";
    }
    close $FH_INFO;

    ## DEBUG
    #~ say "INFO SORTING PRE";
    #~ print Dump @info;
    #my @sorted_info = sort { $a <=> $b } @info;
    #~ say "INFO SORTING POST";
    #~ print Dump @sorted_info;

    my $start = (shift @info);
    my $end   = pop @info;
    ## 2016.12.15a
    ## hack: thresholds may be too strict in some cases???
    ## make sure that at least one position bordering the site is included
    if ($start > $offset - 1)
    {
        $start = $offset - 1;
    }
    ## 2016.12.15a
    ## BUG??
    ## it should be one base after the GT/AG/ATG etc 2fix???
    if ($end < $offset + 1)
    {
        $end = $offset + 1;
    }

    #~ if ($start < 1)
    #~ {
    #~ ## 2016.12.14a
    #~ ## strange check for positions. this should never happen...
    #~ $start = 1;
    #~ }

    print {*STDERR} "\nbit_score_graph out: $start, $end\n";
    ### [<now>] Finished bit_score_graph
    return ($start, $end);
}    #end BitScoreGraph

sub get_pre_matrix ($kmers_tbl, $order)
{
    ### [<now>] Running get_pre_matrix...
    my @orders  = (qw(hmm_0 hmm_2 hmm_3 hmm_4 hmm_5 hmm_6 hmm_7 hmm_8));
    my $ordname = $orders[$order];
    my $FH_KMERS;
    open($FH_KMERS, "<", "$kmers_tbl") or croak "Failed here";
    my $first_line = <$FH_KMERS>;
    my ($seq_name, $seq) = split $first_line;
    my $seq_len = length($seq);
    close $FH_KMERS;

    my $frequency_fn = $work_dir . basename($kmers_tbl) . ".freq";
    my $my_command   = ("./bin/frequency.pypy 1 $kmers_tbl  >  $frequency_fn");
    run($my_command);    ### [<now>] Running frequency.pypy at <file>[<line>]...
    ### [<now>] Finished get_pre_matrix...
}

## GETKMATRIX FUNCTION (Splice sites and Start ATG codon PWMs)
sub get_K_matrix ($true_kmers_tbl, $backgrnd_kmers_tbl, $order, $offset,
                  $matrix_type)
{

    ### [<now>] Running get_K_matrix at <file>[<line>]......
    ## my $original_offset = $offset;
    my @profile_array = ();
    ## $temp_infolog;
    my @orders  = (qw(hmm_0 hmm_2 hmm_3 hmm_4 hmm_5 hmm_6 hmm_7 hmm_8));
    my $ordname = $orders[$order];
    my $sort    = "sort -n";
    ## 2016.12.13a {
    ## BUG do not relay on  jacknife/branch etc. use some: type ??
    #was my our?
    ## my ($true_kmers_tbl, $backgrnd_kmers_tbl, $order, $offset, $matrix_type) =
    ##  @_;

    my ($donor, $acceptor_mtype, $ATGx) = (0, 0, 0);
    ## BUG temp fix
    #~ my $run_jacknife_flag = 0;
    if ($matrix_type eq 'donor')
    {
        $donor = 1;
    }
    elsif ($matrix_type eq 'acceptor')
    {
        $acceptor_mtype = 1;
    }
    elsif ($matrix_type eq 'ATGx')
    {
        $ATGx = 1;    ### [<now>] $matrix_type ATGx...
    }
    else
    {
        croak "sub_get_K_matrix: wrong matrix type\n";
    }

    #     $matrix_type = 'donor';
    # }
    #$branch, $start, $end, $run_jacknife_flag

    # my $matrix_type = "";
    # if ($donor) {
    #     $matrix_type = 'donor';
    # }
    # if ($acceptor_mtype) {
    #     $matrix_type = 'acceptor';
    # }

    # if ($ATGx) {
    #     $matrix_type = 'ATGx';
    # }

    #~ if ($branch) {
    #~ $matrix_type = 'branch';

    #~ }

    if ($order > 1)
    {
        $sort = "sort -k1,1n -k2,2";
    }

    #    my @info = ($offset-1,$offset+1);
    my $profile_len = 0;
    ## 2016.12.15a
    ## my $info_thresh = "";    #bits

    ## BUG?
    my $true_seq_name = $true_kmers_tbl;
    $true_seq_name =~ s/\.tbl$//;
    my $backgrnd_seq_name = $backgrnd_kmers_tbl;
    $backgrnd_seq_name =~ s/\.tbl$//;
    ### true_seq_name: $true_seq_name
    ### backgrnd_seq_name: $backgrnd_seq_name

    ## Open true sequences
    #    print {*STDERR} "$true_kmers_tbl (true)\n";
    open(my $FH_TRUE_SEQ, "<", "$true_kmers_tbl") or croak "Failed here";
    $_ = <$FH_TRUE_SEQ>;
    my @columns_t = split;
    my $len       = length($columns_t[1]);
    close $FH_TRUE_SEQ;

    ## Open false (background???) sequences
    #    print {*STDERR} "$backgrnd_kmers_tbl (false)\n";
    #2017.07.07 skipping: precompute base freqs DK
    #~ open(my $FH_BACKGRND_SEQ, "<", "$backgrnd_kmers_tbl")
    #~ or croak "Couldn't open $backgrnd_kmers_tbl: $OS_ERROR \n";
    #~ $_ = <$FH_BACKGRND_SEQ>;
    #~ my @columns_f = split;
    #~ my $len2      = length($columns_f[1]);
    #~ close $FH_BACKGRND_SEQ;
    my $len2 = 60;

    #    die "$len != $len2\n" if $len != $len2;
    my $true_seq_freq_fn = $work_dir . basename($true_seq_name) . ".freq";
    my $backgrnd_seq_freq_fn =
      $work_dir . basename($backgrnd_seq_name) . ".freq";

    #my $subtracted_true_false_freq_fn = $work_dir . basename($true_seq_name) . "_" . basename($backgrnd_seq_name).freq_subtr";
    my $my_freq_subtract_fn =
        $work_dir
      . basename($true_seq_name) . "_"
      . basename($backgrnd_seq_name)
      . ".information";

    run("./bin/frequency.pypy 1 $true_kmers_tbl  >  $true_seq_freq_fn");

    ### [<now>] not running background freq again at <file>[<line>]...
    ## run("./bin/frequency.pypy 1 $backgrnd_kmers_tbl >  $backgrnd_seq_freq_fn");

    my $my_command_A =
      "./bin/information.pypy  $true_seq_freq_fn $backgrnd_seq_freq_fn ";

    #~ my $my_command_B =
    #~ "| gawk 'NF==2 && \$1<=$my_freq_field_limit_1 && \$1>=$my_freq_field_limit_2'";

    my $my_command_C = " > $my_freq_subtract_fn ";

    #my $my_command = $my_command_A . $my_command_B;
    my $my_command = $my_command_A . $my_command_C;

    say "\n $my_command \n";
    run($my_command);
    ##$temp_infolog = $my_freq_subtract_fn;

    ## True_False req_matrix_fn
    my $my_true_freq_matrix_fn =
      $work_dir . basename($true_seq_name) . "_" . "$ordname.matrix";
    ## False_True req_matrix_fn
    my $my_backgrnd_freq_matrix_fn =
      $work_dir . basename($backgrnd_seq_name) . "_" . "$ordname.matrix";

    ## True logratio req_matrix_fn
    my $my_T_generic_logratio_freq_matrix_fn =
      $work_dir . basename($true_seq_name) . "$matrix_type.log.$ordname.matrix";

    if (!$order)
    {
        $my_command =
          ##"gawk -f ./bin/logratio_zero_order.awk $backgrnd_seq_freq_fn $true_seq_freq_fn > $my_T_generic_logratio_freq_matrix_fn";
          "./bin/logratio_zero_order.pypy $backgrnd_seq_freq_fn $true_seq_freq_fn > $my_T_generic_logratio_freq_matrix_fn";
        say "\n $my_command \n";
        run($my_command);

    }
    else
    {
        $my_command =
          ##"gawk -f ./bin/Getkmatrix.awk $order $len $true_kmers_tbl | $sort > $my_true_freq_matrix_fn";
          "./bin/get_k_matrix.pypy $order $len $true_kmers_tbl | $sort > $my_true_freq_matrix_fn";
        say "\n $my_command \n";
        run($my_command);

        #~ run(
        #~ " gawk -f ./bin/Getkmatrix.awk $order $len $true_kmers_tbl | $sort > $true_seq_name.$ordname-matrix"
        #~ );

        $my_command =
          ##"gawk -f ./bin/Getkmatrix.awk $order $len2 $backgrnd_kmers_tbl | $sort > $my_backgrnd_freq_matrix_fn ";
          "./bin/get_k_matrix.pypy $order $len2 $backgrnd_kmers_tbl | $sort > $my_backgrnd_freq_matrix_fn ";
        say "\n $my_command \n";
        run($my_command);

        $my_command =
          ## "gawk -f ./bin/logratio_kmatrix.awk $my_backgrnd_freq_matrix_fn $my_true_freq_matrix_fn > $my_T_generic_logratio_freq_matrix_fn ";
          "./bin/logratio_kmatrix.py $my_backgrnd_freq_matrix_fn $my_true_freq_matrix_fn > $my_T_generic_logratio_freq_matrix_fn ";

        say "\n $my_command \n";
        run($my_command);

    }

    #need to check output and then go on
    ## draw bit score bar graph function (nested, local)

    #~ local *BitScoreGraph = sub {

    #~ my ( $info_output, $info_thresh, $offset ) = @_;
    #~ my @info = ( $offset - 1, $offset + 1 );
    #~ my @fields;
    #~ open( my $FH_INFO, "<", "$info_output" ) or croak "Failed here";
    #~ while (<$FH_INFO>) {
    #~ next if m/^#/;
    #~ last if m/^\s/;
    #~ last if m/^[^\d]/;
    #~ chomp;
    #~ @fields = split;
    #~ printf {*STDERR} "%2s %2.2f %s",
    #~ ##
    #~ ( $fields[0], $fields[1], "=" x int( $fields[1] * 30 ) );
    #~ if ( $fields[1] > $info_thresh ) {
    #~ push( @info, $fields[0] );
    #~ }
    #~ print {*STDERR} "\n";
    #~ }
    #~ close $FH_INFO;
    #~ print {*STDERR} "\n BitScoreGraph \n";

    #~ my @sorted_info = sort numerically @info;
    #~ my $start      = ( shift @sorted_info );
    #~ if ( $start < 1 ) {
    #~ $start = 1;
    #~ }
    #~ my $end = pop @sorted_info;

    #~ return ( $start, $end );
    #~ };    #end BitScoreGraph

    ## TODO simplify

    #~ $run_jacknife_flag = 0;

    #~ if ( !$run_jacknife_flag ) {
    #~ print {*STDERR} "Information content profile\n";
    #~ }

    #~ my %my_info_thresholds = (
    #~ donor    => 0.15,
    #~ acceptor => 0.04,
    #~ ATGx     => 0.15,
    #~ branch   => 0.30,
    #~ );
    print {*STDERR}
      "\n L2917: $matrix_type, $my_info_thresholds{$matrix_type}\n";

    my $my_info_thresh = $my_info_thresholds{$matrix_type};

    #say "\n matrix sub: $matrix_type, $my_info_thresh \n";
    my ($start, $end) =
      bit_score_graph($my_freq_subtract_fn, $my_info_thresh, $offset);
    print {*STDERR} "\n L2921 got: $start, $end using $my_freq_subtract_fn \n";

    $offset = $offset - $order;
    $end    = $end - $order;

    #    if ( $start < 1 ) {
    #        $start = 1;
    #    }
    ## BUG changing offset again!!!
    $offset = $offset - $start + 1;
    print {*STDERR}
      "end:$end offset:$offset start:$start  donor:$donor acceptor_mtype:$acceptor_mtype ATG: $ATGx\n";
    print {*STDERR}
      "new offset: $offset\nnew start: $start\nnew order: $order\n";

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
    my $exec_A1 = "./bin/submatrix.py ";
    ##my $exec_B1 = "./bin/preparedimatrixdonor4parameter.awk ";
    my $exec_B1 = "./bin/matrix_4_parameter.py ";

    if ($order >= 1 && $donor)
    {
        ### [<now>] Running donor order_oneplus...
        #my $pre_offset = $offset + $my_dimatrics_dict{'donor'};

        my $pre_offset  = $offset + 2;
        my $new_offset  = $offset + 3;
        my $post_offset = $offset + 4;

        #my_True_dimatrixdonor_4param_fn
        #my $exec_A1 = "gawk -f ./bin/submatrix.awk ";

        my $my_command_A =
          "$exec_A1 $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_lograt_summatrix_fn";

        #~ my $my_command_B =
        #~ "$exec_B1 $pre_offset $new_offset $post_offset $my_T_generic_lograt_summatrix_fn > $my_T_generic_matrix_4param_fn";

        my $my_command_B =
          "$exec_B1 $my_T_generic_lograt_summatrix_fn donor $pre_offset $new_offset $post_offset  > $my_T_generic_matrix_4param_fn";

        print {*STDERR} "$my_command_A \n";
        run($my_command_A);
        print {*STDERR} "$my_command_B \n";
        run($my_command_B);

        # print {*STDERR} "submatrix.awk $start $end $true_seq_name-log.$ordname-matrix/$true_seq_name-log-info.$ordname-matrix";

    }
    ## DONOR DIMATRIX END
    ## ACCEPTOR DIMATRIX START
    elsif ($order >= 1 && $acceptor_mtype)
    {
        ### [<now>] Running acceptor order_one_plus...
        my $pre_offset  = $offset - 1;
        my $new_offset  = $offset;
        my $post_offset = $offset + 1;

        ## my $my_command_A =
        ## "gawk -f ./bin/submatrix.awk $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_lograt_summatrix_fn";
        my $my_command_A =
          "$exec_A1 $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_lograt_summatrix_fn";

        print {*STDERR} "my_command_A: \n$my_command_A \n";
        run($my_command_A);

        #~ my $my_command_B =
        #~ "./bin/preparedimatrixacceptor4parameter.awk $pre_offset $new_offset $post_offset $my_T_generic_lograt_summatrix_fn > $my_T_generic_matrix_4param_fn";

        my $my_command_B =
          "$exec_B1 $my_T_generic_lograt_summatrix_fn acceptor $pre_offset $new_offset $post_offset  > $my_T_generic_matrix_4param_fn";

        print {*STDERR} "my_command_B: \n$my_command_B \n";
        run($my_command_B);

        #~ run(
        #~ "./bin/preparedimatrixacceptor4parameter.awk $pre_offset $new_offset $post_offset $my_Tlograt_summatrix_fn > $my_T_dimatrixdonor_4param_fn"
        #~ );

        #~ #      print {*STDERR} "submatrix.awk $start $end $true_seq_name-log.$ordname-matrix/$my_True_dimatrixdonor_4param_fn";

    }
    ## ACCEPTOR DIMATRIX END

    ## ATG DIMATRIX START
    elsif ($order >= 2 && $ATGx)
    {
        ### [<now>] Running ATGx order_two...
        my $pre_offset  = $offset - 2;
        my $new_offset  = $offset - 1;
        my $post_offset = $offset;
        $my_command =
          "./bin/submatrix.py $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_lograt_summatrix_fn";
        ##run("gawk -f ./bin/submatrix.awk $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_lograt_summatrix_fn"
        ##);
        say $my_command;
        run($my_command);

        run(" ./bin/preparetrimatrixstart4parameter.awk $pre_offset $new_offset $post_offset $my_T_generic_lograt_summatrix_fn > $my_T_generic_matrix_4param_fn"
           );
    }

    ## ATG DIMATRIX END

    ## ALL REMAINING CASES START
    else
    {    ### [<now>] Running REMAINING CASES...

        # print {*STDERR} "$path/submatrix_order0.awk $start $end $true_seq_name-log.$ordname-matrix\n";
        $my_command =
          "./bin/submatrix_order0.py $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_matrix_4param_fn";
        say $my_command;
        run($my_command);
        ## run(
        ##" gawk -f ./bin/submatrix_order0.awk $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_matrix_4param_fn"
        ##        );

    }
    ## ALL REMAINING CASES END

    ## CREATE DATA STRUCTURE CONTAINING MATRIX OF INTEREST

    open(my $FH_PROF, "<", "$my_T_generic_matrix_4param_fn")
      or croak "Failed here";
    while (<$FH_PROF>)
    {
        next if m/^#/;
        last if m/^\s/;
        last if m/^[^\d]/;
        chomp;
        my @columns_e = split;
        push @profile_array, \@columns_e;
    }
    close $FH_PROF;

    $profile_len = $end - $start + 1;
    print {*STDERR} "??profile??  length: $end - $start / $profile_len \n";
    return (\@profile_array, $profile_len, $offset, $start, $end);

    #unlink $temp_infolog;
    #unlink "$my_T_dimatrixdonor_4param_fn";
    #unlink "$my_Tlograt_summatrix_fn";

}    ##END GETKMATRIX FUNCTION (Splice site an Start codon PWMs)

## sub numerically { $a <=> $b }

sub fold4fasta ($seq)
{
    ## 2016.12.13a {
    ## my $seq        = shift;
    my $folded_seq = "";

    #my $column_limit = 60;
    Readonly::Scalar my $column_limit => 60;

    #my $s = "";
    my $position_in_seq = 0;
    my $seq_len         = length($seq);

    while ($position_in_seq < $seq_len)
    {
        my $out_seq = substr($seq, $position_in_seq, $column_limit);
        $folded_seq = $folded_seq . $out_seq . "\n";
        $position_in_seq += $column_limit;
    }
    return $folded_seq;
}

#Optimize parameter file
sub parameter_optimize ($gp_fasta, $gp_gff_fn, $new_param_fn, %profile_params)
{

    ## 2016.12.13a
    ## my ($gp_fasta, $gp_gff_fn, $new_param_fn, %profile_params) = @_;
    ## TODO these 4 vars were not used in that sub
    #        $branch_switch,
    #        $branch_profile_len,
    #        $fxdbraoffset,
    #        $branch_matrix,

    my $ExWeightParam_ini   = $profile_params{'ExWeightParam_ini'};
    my $ExWeightParam_step  = $profile_params{'ExWeightParam_step'};
    my $ExWeightParam_limit = $profile_params{'ExWeightParam_limit'};
    ## EXON/OLIGO FACTOR PARAMETER
    my $OligoWeight_ini   = $profile_params{'OligoWeight_ini'};
    my $OligoWeight_step  = $profile_params{'OligoWeight_step'};
    my $OligoWeight_limit = $profile_params{'OligoWeight_limit'};
    ## Minimum Branch Profile Distance
    ## 20161201 unused
    ## my $i_BranchProfDist = $profile_params{'i_BranchProfDist'};
    ## my $d_BranchProfDist = $profile_params{'d_BranchProfDist'};
    ## my $f_BranchProfDist = $profile_params{'f_BranchProfDist'};

    ## my $iAccCtx = $profile_params{'iAccCtx'};
    ## my $dAccCtx = $profile_params{'dAccCtx'};
    ## my $fAccCtx = $profile_params{'fAccCtx'};

    my @evaluation_total    = ();
    my $myExWeightParam_tmp = $ExWeightParam_ini;
    my $myOlWeight_tmp      = $OligoWeight_ini;
    ##unused 20161201
    ##my $i_BranchProfDistini          = $i_BranchProfDist;
    ##my $iAccCtxini       = $iAccCtx;

    #my $FH_SOUT;

    open(my $FH_SOUT, ">", "$work_dir/$species##.OptimizeParameter.log")
      or croak "Failed here";

    print {*STDERR}
      "\neWF range : $ExWeightParam_ini to $ExWeightParam_limit\noWF range : $OligoWeight_ini to $OligoWeight_limit\n\n";
    print {$FH_SOUT}
      "\neWF range : $ExWeightParam_ini to $ExWeightParam_limit\noWF range : $OligoWeight_ini to $OligoWeight_limit\n\n";
    close $FH_SOUT;

    for ($myExWeightParam_tmp = $ExWeightParam_ini ;
         $myExWeightParam_tmp <= $ExWeightParam_limit ;
         $myExWeightParam_tmp += $ExWeightParam_step)
    {    #for_#1
        print {*STDERR} "eWF: $myExWeightParam_tmp\noWF: ";

        for ($myOlWeight_tmp = $OligoWeight_ini ;
             $myOlWeight_tmp <= $OligoWeight_limit ;
             $myOlWeight_tmp += $OligoWeight_step)
        {    #for_#2
            print {*STDERR} "$myOlWeight_tmp";
            my $param = Geneid::Param->new();
            $param->readParam("$new_param_fn");

            ##20170222 for (my $i = 0 ; $i < $param->numIsocores ; $i++)
            for my $i (0 .. $param->numIsocores - 1)
            {    #for_#3
                if (
                    !defined @{$param->isocores}[$i]->Exon_weights(
                                 [
                                  $myExWeightParam_tmp, $myExWeightParam_tmp,
                                  $myExWeightParam_tmp, $myExWeightParam_tmp
                                 ]
                    )
                   )
                {
                    croak "error in setting exon weights L3163\n";
                }

                if (
                    !defined @{$param->isocores}[$i]->Exon_factor(
                                           [
                                            $myOlWeight_tmp, $myOlWeight_tmp,
                                            $myOlWeight_tmp, $myOlWeight_tmp
                                           ]
                    )
                   )
                {
                    croak "error in setting exon weights\n";
                }

                #~ #   if (!defined @{$param->isocores}[$i]->Exon_factor([0.4,$IoWF,$IoWF,0.4])) {

                if (
                    !defined @{$param->isocores}[$i]->Site_factor(
                                                        [
                                                         1 - $myOlWeight_tmp,
                                                         1 - $myOlWeight_tmp,
                                                         1 - $myOlWeight_tmp,
                                                         1 - $myOlWeight_tmp
                                                        ]
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
            #
            #~ my $FH_GENEID    = File::Temp->new();
            #~ my $fname_geneid = $FH_GENEID->filename;
            my $FH_GENEID;
            my $fname_geneid;
            ($FH_GENEID, $fname_geneid) =
              tempfile(DIR => $geneid_dir, SUFFIX => '.geneid');
            print "\ntemp geneid file: $fname_geneid \n";
            my $my_command =
              "./bin/geneid -GP $temp_geneid_param $gp_fasta > $fname_geneid";
            run($my_command);
            my $temp_geneid_pred_gff_fn =
              "$tmp_dir/Predictions." . basename($new_param_fn) . ".gff";
            $my_command =
              "cat $fname_geneid | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' > $temp_geneid_pred_gff_fn ";
            run($my_command);

            #` ./bin/geneid -GP ${newparam}.temp $gp_fasta | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' > $tmp_dir/Predictions.${newparam}.gff`;
            ## BUG very complex comand line

            my $temp_evalout_A_fn =
                "$tmp_dir/Predictions."
              . basename($new_param_fn)
              . ".temp_evalout_A";
            $my_command =
              "./bin/evaluation -sta $temp_geneid_pred_gff_fn $gp_gff_fn  > $temp_evalout_A_fn";
            print "\n$my_command\n";
            run($my_command);

            my $temp_evalout_B_fn =
                "$tmp_dir/Predictions."
              . basename($new_param_fn)
              . ".temp_evalout_B";
            $my_command = "tail -2 $temp_evalout_A_fn | head -1 |
            gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $myOlWeight_tmp, $myExWeightParam_tmp, \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}' > $temp_evalout_B_fn ";
            print "\n$my_command\n";
            run($my_command);

            my @evaluation_output;
            open(my $FH_IN, "<", "$temp_evalout_B_fn") or croak "Failed here";
            while (<$FH_IN>)
            {
                #                @evaluation_output = split " ";
                @evaluation_output = split;

                #@evaluation_output = split " ", $_;
            }
            close $FH_IN;

            #

            #~ my @evaluation_output = split " ",

            #~ #` ./bin/evaluation -sta $temp_geneid_pred_gff_fn $gp_gff_fn | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $IoWF, $I_ExWeight_F, \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}' `;

            push(@evaluation_total, \@evaluation_output);

        }    #end for_#2
        $myOlWeight_tmp = $OligoWeight_ini;
        print {*STDERR} "\n";

    }    #end for_#1

    return \@evaluation_total;

}
## end sub optimize parameter file

sub get_opt_paramfile ( $eval_array )
{
    ### [<now>] Running get_opt_paramfile...
    ##,   $branch_switch, $branch_profile_len,
    ##$fxdbraoffset, $branch_matrix
    ##)

    ## 2016.12.13a
    #~ my (
    #~ $eval_array,   $branch_switch, $branch_profile_len,
    #~ $fxdbraoffset, $branch_matrix
    #~ )
    #~ = @_;
    my @sorted_eval        = ();
    my @evaluation_init    = ();
    my $best_OlWeight      = 0;
    my $best_ExWeightParam = 0;
    ## 2016.12.15 unused
    ## my $best_Min           = 0;
    ## my $best_Acc           = 0;
    ## $FH_SOUT;
    open(my $FH_SOUT, ">", "$work_dir/$species.get_opt_paramfile.log")
      or croak "Failed here";

    ## DEBUG

    #~ print {*STDERR} "\n input to OptimizeParameter sub
    #~ \n1:"
    #~ . $eval_array . "\n2 "
    #~ . $branch_switch . "\n3 "
    #~ . $branch_profile_len . "\n4 "
    #~ . $fxdbraoffset . "\n5 "
    #~ . $branch_matrix . "\n";

    ##if (!$branch_switch)
    ##if (1 > 0)
    ## {

    ## DEBUG
    #~ say "EVAL_ARRAY";
    #~ print Dump $eval_array;
    ## BUG ???
    @sorted_eval = sort sorteval @{$eval_array};

    #~ say "EVAL_ARRAY_SORTED";
    #~ print Dump @sorted_eval;

    $best_OlWeight      = $sorted_eval[0][0];    #0.2
    $best_ExWeightParam = $sorted_eval[0][1];    #-3.5

    print {*STDERR} "\nBest performance obtained using I_OlWeight_F: "
      . $sorted_eval[0][0]
      . " and myExWeightParam_tmp: "
      . $sorted_eval[0][1] . "\n";
    print {$FH_SOUT} "\nBest parameter file performance obtained using oWF: "
      . $sorted_eval[0][0]
      . " and eWF: "
      . $sorted_eval[0][1] . "\n";

    #INITIALIZE ARRAY WITH EVALUATION PARAMETERS
    @evaluation_init =
      (qw(oWF eWF SN SP CC SNe SPe SNSP SNg SPg SNSPg raME raWE));

    print {*STDERR}
      "Sorted performance results (Three best performance estimates) for different values of oWF and eWF:\n"
      . join("\t", @evaluation_init), "\n";
    print {$FH_SOUT}
      "Sorted performance results (best to worst) for different values of oWF and eWF:\n\n"
      . join("\t", @evaluation_init), "\n";

    foreach my $eval_ref (@sorted_eval)
    {

        print {$FH_SOUT} join("\t", @{$eval_ref}), "\n";

    }

    ## FOUR BEST PERFORMANCE RESULTS AFTER OPTIMIZATION
    ## BUG??? 0, 1, 3 is just 3 results ??
    ##for (my $i = 0 ; $i <= 2 ; $i++)
    for my $ij (0 .. 2)
    {
        print {*STDERR} join("\t", @{$sorted_eval[$ij]}), "\n";
    }
    ##

    ## BUILD BEST PERFORMANCE (OPTIMIZED) PARAMETER FILE (USE VALUES OBTAINED ABOVE)

    my $param = Geneid::Param->new();
    $param->readParam("$new_param_fn");

    for (my $kk = 0 ; $kk < $param->numIsocores ; $kk++)
      ##for my $kk (0 .. $param->numIsocores )
    {
        if (
            !defined @{$param->isocores}[$kk]->Exon_weights(
                                   [
                                    $best_ExWeightParam, $best_ExWeightParam,
                                    $best_ExWeightParam, $best_ExWeightParam
                                   ]
            )
           )
        {
            croak "error in setting exon weights\n";
        }
        if (
            !defined @{$param->isocores}[$kk]->Exon_factor(
                                             [
                                              $best_OlWeight, $best_OlWeight,
                                              $best_OlWeight, $best_OlWeight
                                             ]
            )
           )
        {
            #     if (!defined @{$param->isocores}[$kk]->Exon_factor([0.4,$best_IoWF,$best_IoWF,0.4])) {
            croak "error in setting exon weights\n";
        }
        if (
            !defined @{$param->isocores}[$kk]->Site_factor(
                                                         [
                                                          1 - $best_OlWeight,
                                                          1 - $best_OlWeight,
                                                          1 - $best_OlWeight,
                                                          1 - $best_OlWeight
                                                         ]
            )
           )
        {
            #     if (!defined @{$param->isocores}[$kk]->Site_factor([0.55,1-$best_IoWF,1-$best_IoWF,0.55])) {
            croak "error in setting exon weights\n";
        }
    }

    #write new parameter file (optimized)
    #    my $optimized_geneid_param_fn =
    #      "$results_dir/$species.geneid.optimized.param";
    my $optimized_geneid_param_fn = "$species.geneid.optimized.param";
    $param->writeParam($optimized_geneid_param_fn);

    #    print {*STDERR}
    #      "\n\nNew optimized parameter file named: $species.geneid.optimized.param \n";
    #    print {$FH_SOUT}
    #      "\n\nNew optimized parameter file named: $species.geneid.optimized.param \n";

    close $FH_SOUT;
    ### [<now>] Finished get_opt_paramfile...
    ### $optimized_geneid_param_fn

    return [$best_ExWeightParam, $best_OlWeight, 0, 0, \@evaluation_init];

}

sub parameter_evaluate ($gp_fasta, $gp_gff_fn, $new_param_fn, $OligoWeight_ini,
                        $ExWeightParam_ini)
{
    ### [<now>] Running  parameter_evaluate...

    my $my_command;
    my $geneid_test_predict_gff_fn =
      "$work_dir/geneid_test_predictions." . basename($new_param_fn) . ".gff";
    ### $geneid_test_predict_gff_fn
    #    print {*STDERR}
    #      "\ngeneid_test_predict_gff_fn: $geneid_test_predict_gff_fn\n";

    my ($FH_GENEID, $fname_geneid) =
      tempfile(DIR => $geneid_dir, SUFFIX => '.eval_par.geneid');
    print "\ntemp geneid file: $fname_geneid \n";
    $my_command = "./bin/geneid -GP $new_param_fn $gp_fasta > $fname_geneid";

    #print {*STDERR} "\n$my_command, not running\n";
    run($my_command);

    #` ./bin/geneid -GP $new_param_fn $gp_fasta
    $my_command =
      "cat $fname_geneid | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' > $geneid_test_predict_gff_fn";
    run($my_command);

    #dk
    #` ./bin/geneid -GP ${newparam}.temp $gp_fasta | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' > $tmp_dir/Predictions.${newparam}.gff`;
    ## BUG very complex comand line

    my $temp_evalout_A_fn =
      "$geneid_dir/Predictions." . basename($new_param_fn) . ".temp_evalout_A";
    $my_command =
      "./bin/evaluation -sta $geneid_test_predict_gff_fn $gp_gff_fn > $temp_evalout_A_fn";
    print "\n$my_command\n";
    ### [<now>] evaluation
    ### $my_command
    run($my_command);

    my $temp_evalout_B_fn =
      "$geneid_dir/Predictions." . basename($new_param_fn) . ".temp_evalout_B";
    $my_command = "tail -2 $temp_evalout_A_fn | head -1 |
            gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $OligoWeight_ini, $ExWeightParam_ini, \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}' > $temp_evalout_B_fn ";
    print "\n$my_command\n";
    run($my_command);

    my @evaluation_test;
    open(my $FH_IN, "<", "$temp_evalout_B_fn") or croak "Failed here";
    while (<$FH_IN>)
    {
        #@evaluation_test = split " ";
        @evaluation_test = split;

        #@evaluation_output = split " ", $_;
    }
    close $FH_IN;

    #dk

    #~ my @evaluation_test = split " ",
    #~ ` ./bin/evaluation -sta $geneid_test_predict_gff_fn $gp_gff_fn | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}' `;
    ### [<now>] Finished sub  parameter_evaluate...
    return \@evaluation_test;

}    # evaluate parameter function

sub calc_stats
{
    ## BUG variable names hard to guess
    my (
        ## 2016.12.15 unused
        ## $species,
        ## $sout,
        $train_introns_filtered_tbl, $train_cds_filtered_tbl,
        $train_filtered_gff,         $train_inframestop_int,
        $eval_inframestop_int,       $train_transcr_used_int,
        $train_noncanon_donors_int,  $train_noncanon_accept_int,
        $train_noncanon_ATGx_int,    $markov_model,
        $total_coding,               $total_noncoding,
        $st_donor,                   $en_donor,
        $st_accept,                  $en_accept,
        $st_ATGx,                    $en_ATGx,

        #$st_branch,
        #$en_branch,
        #$branch_switch,
        #$use_allseqs_flag
       ) = @_;
    my $my_command = "";

    ### [<now>] Running calc_stats...
    ## OBTAIN GENE MODEL SET STATISTICS
    ## Open gene model object

    $param->geneModel(Geneid::GeneModel->new());
    $param->geneModel->useDefault;
    #
    ##my $FH_SOUT;
    open(my $FH_SOUT, ">", "$work_dir/test.WriteStatsFileLog.txt")
      or croak "Failed here";

    ## my $avg_intron = "";
    my $sd_intron = "";

    #my @intronlength = ` gawk '{print length(\$2)}' $introns_clean_tbl_fn `;

    #print {*STDERR} "INTRON: $mean, $st\n";
    $my_command =
      "gawk '{print length(\$2)}' $train_introns_filtered_tbl | sort -n";
    my @introns_len_list = capture($my_command);
    my ($mean, $st) = calc_average(\@introns_len_list);
    print {*STDERR} "INTRONS mean, ST: $mean, $st\n";

    #~ ## BUG wrong command?

    #~ #$my_command = "gawk '{print length(\$2)}' $introns_clean_tbl_fn | sort -n";
    #~ $my_command = "sort -k2,2n $introns_clean_tbl_fn";
    #~ my @intronlist = capture($my_command);

    ## my $total_introns     = scalar(@introns_len_list);
    my @intron_lenX_array = ();
    my $intron_len        = "";
    for (my $i = 0 ; $i <= scalar(@introns_len_list) - 1 ; $i++)
    {

        $intron_len = $introns_len_list[$i];
        chomp $intron_len;
        push(@intron_lenX_array, $intron_len);
    }

    my @slice1 = @intron_lenX_array[0 .. 5];
    my @slice2 =
      @intron_lenX_array[(scalar(@intron_lenX_array) - 5)
      .. (scalar(@intron_lenX_array) - 1)];
    ##@intron_lenX_array[ ( scalar(@intron_lenX_array) - 5 ) .. ( scalar(@intron_lenX_array) - 1 ) ]
    ## BUG 2857
    my $intron_short_int = $introns_len_list[0] - $introns_len_list[0] * (0.25);
    chomp $intron_short_int;
    if ($intron_short_int > 40) { $intron_short_int = 40; }

    #my $intron_long_int =  $intronlist[$total_introns - 1];
    my $intron_long_int =
      $mean + ($st * 3) > 100_000 ? 100_000 : $mean + ($st * 3);
    chomp $intron_long_int;

    my $intergenic_min = 200;
    my $intergenic_max = 'Infinity';

    ## use shortest and longest intron lengths in gene model of parameter file
    $param->geneModel->intronRange($intron_short_int, $intron_long_int);
    $param->geneModel->intergenicRange($intergenic_min, $intergenic_max);
    #

    $my_command =
      "gawk '{print gsub(/[GC]/,\".\",\$2)/length(\$2)}' $train_cds_filtered_tbl";
    my @cds_GC_content = capture($my_command);

    #print {*STDERR} "@cds_GC_content\n";
    my ($mean_GC_CDS, $stdev_GC_CDS) = calc_average(\@cds_GC_content);

    #print {*STDERR} "CDS: $meangc $stgc $introns_clean_tbl_fn\n";
    $my_command =
      "gawk '{print gsub(/[GC]/,\".\",\$2)/length(\$2)}' $train_introns_filtered_tbl ";
    my @intron_GC_content = capture($my_command);

    #print {*STDERR} "@intron_GC_content\n";
    my ($mean_GC_intron, $stdev_GC_intron) = calc_average(\@intron_GC_content);

    #print {*STDERR} "intron: $meangci $stgci\n";
    #BUG?
    #my $total_exons = ` gawk '{print \$9}' $train_filtered_gff | wc -l | gawk '{print \$1}' `;
    $my_command = "gawk '{print \$9}' $train_filtered_gff | wc -l ";
    my $total_exons = capture($my_command);

    #my $total_exons = ` gawk '{print \$9}' $train_filtered_gff | wc -l `;

    chomp $total_exons;
    $total_exons = int($total_exons);
    my @exons_per_gene;
    $my_command =
      "gawk '{print \$9}' $train_filtered_gff | sort | uniq -c | gawk '{print \$1}'";
    @exons_per_gene = capture($my_command);

    #@exons_per_gene =
    #  ` gawk '{print \$9}' $train_filtered_gff | sort | uniq -c | gawk '{print \$1}' `;

    my ($avg_ex, $stdevX_ex) = calc_average(\@exons_per_gene);

    $my_command =
      "egrep -v 'Single' $train_filtered_gff | gawk '{len=\$5-\$4;print len}' - | sort ";
    my @exons_lenghts_list = capture($my_command);

    #    my @exons_lenghts_list =
    #      ` egrep -v 'Single' $train_filtered_gff | gawk '{len=\$5-\$4;print len}' - | sort `;

    my ($avg_exon_len, $stdev_exonX_len) = calc_average(\@exons_lenghts_list);

    $my_command = "egrep -c '(Single)' $train_filtered_gff";
    my $single_exon_genes = capture($my_command);

    #my $single_exon_genes = `egrep -c '(Single)' $train_filtered_gff `;
    chomp $single_exon_genes;

    #print {$FH_SOUT} "GENE MODEL STATISTICS FOR $species\n\n";

    print {$FH_SOUT}
      "\nA subset of $train_transcr_num sequences (randomly chosen from the $transcr_all_number gene models) was used for training\n\n";
    ## my $transcr_all_number;

    #~ if (!$use_allseqs_flag)
    #~ {
    print {$FH_SOUT}
      "The user has selected to use $train_transcr_used_int gene models (80 % of total) for training and to set aside BUG  was xxgffseqseval annotations (20 % of total) for evaluation\n\n";

    #~ }
    #~ else
    #~ {
    #~ print {$FH_SOUT}
    #~ "$transcr_all_number gene models were used for both training and evaluation\n\n";
    #~ }

    #~ if (!$use_allseqs_flag)
    #~ {
    print {$FH_SOUT}
      "$train_inframestop_int of the gene models translate into proteins with in-frame stops within the training set and $eval_inframestop_int in the evaluation set (seqs removed).\n\n";

    #~ }

    #~ else
    #~ {
    #~ print {$FH_SOUT}
    #~ "$inframe_X of the gene models translate into proteins with in-frame stops within the training set.\n\n";
    #~ }
    print {$FH_SOUT}
      "There are $train_noncanon_donors_int non-canonical donors as part of the training set\n\n";
    print {$FH_SOUT}
      "There are $train_noncanon_accept_int non-canonical acceptors as part of the training set\n\n";
    print {$FH_SOUT}
      "There are $train_noncanon_ATGx_int non-canonical start sites as part of the training set\n\n";
    print {$FH_SOUT}
      "These gene models correspond to $total_coding coding bases and $total_noncoding non-coding bases\n\n";
    print {$FH_SOUT}
      "Deriving a markov model for the coding potential of order $markov_model\n\n";
    print {$FH_SOUT}
      "The intronic sequences extracted from the gene models have an average length of $mean, with $st of SD\n";
    print {$FH_SOUT}
      "Geneid can predict gene models having introns with a minimum length of $intron_short_int nucleotides and a maximum of $intron_long_int bases (boundaries used in gene model) \n\n";
    print {$FH_SOUT}
      "The minimum (user selected) intergenic distance was set to $intergenic_min nucleotides whereas the maximum was set to $intergenic_max (boundaries used in gene model) \n\n";
    print {$FH_SOUT}
      "The GC content of the exonic (XXX bug ??? XXX) and intronic sequences is $mean_GC_CDS (SD $stdev_GC_CDS) and $mean_GC_intron (SD $stdev_GC_intron) respectively \n\n";
    print {$FH_SOUT}
      "The gene models used for training contain $total_exons exons \n\n";
    print {$FH_SOUT}
      "The gene models average $avg_ex exons per gene (SD $stdevX_ex)\n\n";
    print {$FH_SOUT}
      "The average length of the exons (non-single) in the training set gene models is $avg_exon_len (SD $stdev_exonX_len)\n\n";

    #~ if (!$use_allseqs_flag)
    #~ {
    print {$FH_SOUT}
      "The training set includes $single_exon_genes single-exon genes (out of $train_transcr_used_int ) gene models\n\n";

    #~ }
    #~ else
    #~ {
    #~ print {$FH_SOUT}
    #~ "The training set includes $single_exon_genes single-exon genes (out of $transcr_all_number) gene models\n\n";
    #~ }
    print {$FH_SOUT} "The donor site profile chosen by the user spans "
      . ($en_donor - $st_donor + 1)
      . " nucleotides: position $st_donor to $en_donor\n";
    print {$FH_SOUT} "The acceptor site profile chosen by the user spans "
      . ($en_accept - $st_accept + 1)
      . " nucleotides: position $st_accept to $en_accept\n";
    print {$FH_SOUT} "The start site profile chosen by the user spans "
      . ($en_ATGx - $st_ATGx + 1)
      . " nucleotides: position $st_ATGx to $en_ATGx\n";
    ## 2016.2.13b
    ## if ($branch_switch)
    ## {
    ##    print {$FH_SOUT} "The branch site profile chosen by the user spans "
    ##      . ($en_branch - $st_branch + 1)
    ##      . " nucleotides: position $st_branch to $en_branch\n";
    ##}
    close $FH_SOUT;
    ### [<now>] Finished calc_stats...
    return ($intron_short_int, $intron_long_int, $intergenic_min,
            $intergenic_max);

}

sub calc_average ($input_numbers)
{
    ### [<now>] Running calc_average...
    my $sum            = 0;
    my $elements_count = 0;
    my ($mean, $stdev);

    foreach my $my_number (@{$input_numbers})
    {
        $sum += $my_number;
        $elements_count++;
    }

    $mean = $sum / $elements_count;
    $mean = sprintf("%.3f", $mean);

    $sum = 0;
    my $tmp_number = 0;
    foreach my $my_number (@{$input_numbers})
    {
        $tmp_number = $my_number - $mean;
        $sum += $tmp_number * $tmp_number;

        #$sum += ( $my_number - $mean ) * ( $my_number - $mean );
    }
    $stdev = sqrt($sum / $elements_count);
    $stdev = sprintf("%.3f", $stdev);
    ### [<now>] Finished calc_average...
    return ($mean, $stdev);
}

sub tbl_2_fasta ($in_tbl_fn, $fa_out_fn)
{
    ### [<now>] Running tbl_2_fasta...

    open(my $FH_IN,   "<", "$in_tbl_fn") or croak "Failed here";
    open(my $FH_FOUT, ">", "$fa_out_fn") or croak "Failed here";
    while (<$FH_IN>)
    {
        chomp;

        #~ my ( $n, $s ) = split( /\s+/, $_ );
        my ($seq_name, $seq) = split;
        my ($seq_position, $seq_len) = (1, length($seq));
        print {$FH_FOUT} ">$seq_name\n";
        while ($seq_position <= $seq_len)
        {
            print {$FH_FOUT} substr($seq, $seq_position - 1, 60) . "\n";
            $seq_position += 60;
        }
    }
    close $FH_IN;
    close $FH_FOUT;
    ### [<now>] Finished tbl_2_fasta...
    return $fa_out_fn;
}

sub tbl_2_single_fastas ($in_tbl_fn, $out_fas_dir)
{
    ### [<now>] Running tbl_2_single_fastas...

    open(my $FH_IN_TBL, "<", "$in_tbl_fn") or croak "Failed here";

    #    print {*STDERR} "##tbl_2_single_fastas $in_tbl_fn\n";
    while (<$FH_IN_TBL>)
    {
        #        my $input_line = $_;
        #        chomp($input_line);

        #my ( $seq_name, $seq ) = split( /\s+/o, $_ );
        #        my ($seq_name, $seq) = split(/\s+/o, $input_line);
        my ($seq_name, $seq) = split;

        #print {*STDERR} "YYY $seq_name \t";
        ##open( FOUT, ">${dir}" . "$n" );
        open(my $FH_FOUT_FASTA, ">", "${out_fas_dir}" . "$seq_name");

        #print {*STDERR} "XXX ${dir}" . "$seq_name \n";
        my ($base_number, $seq_len) = (1, length($seq));
        print {$FH_FOUT_FASTA} ">$seq_name\n";
        print {*STDERR} "#";
        while ($base_number <= $seq_len)
        {
            print {$FH_FOUT_FASTA} substr($seq, $base_number - 1, 60) . "\n";
            $base_number += 60;
        }
        close $FH_FOUT_FASTA;
    }

    close $FH_IN_TBL;
    ### [<now>] Finished tbl_2_single_fastas...
    #close $FH_FOUT_FASTA;
    return 1;
}

sub fasta_2_tbl ($in_fa, $out_tbl_fn)
{
    ### [<now>] Running fasta_2_tbl...

    open(my $FH_IN,   "<", "$in_fa");
    open(my $FH_TOUT, ">", "$out_tbl_fn");

    #print {*STDERR} "$in_fa INSIDE LOOP\n\n";
    my $count = 0;
    while (<$FH_IN>)
    {
        chomp;
        $_ =~ s/\|//;
        if ($_ =~ /\>(\S+)/)
        {
            if ($count > 0)
            {
                print {$FH_TOUT} "\n";
            }
            print {$FH_TOUT} $1 . "\t";
            $count++;
        }
        else
        {
            print {$FH_TOUT} $_;
        }
    }
    print {$FH_TOUT} "\n";

    close $FH_IN;
    close $FH_TOUT;
    ### [<now>] Finished fasta_2_tbl...
    return $out_tbl_fn;
}

#~ sub FastaToTbl {

#~ my ( $in_fa_fn, $out_tbl_fn, $flag ) = @_;
#~ say "\n$in_fa_fn, $out_tbl_fn, $flag\n";
#~ #fix CAPS & Ns while converting
#~ #

#~ open( my $FH_IN,   "<", "$in_fa_fn" ) ;
#~ open( my $FH_TOUT, ">", "$out_tbl_fn" );

#~ #print {*STDERR} "$fa INSIDE LOOP\n\n";
#~ my $count    = 0;
#~ my $sequence = "";
#~ while (<$FH_IN>) {
#~ chomp;
#~ $_ =~ s/\|//;
#~ if ( $_ =~ /\>(\S+)/ ) {
#~ print $FH_TOUT "\n" if $count > 0;
#~ print $FH_TOUT $1 . "\t";
#~ $count++;
#~ }
#~ else {
#~ if ($flag eq "genome"){
#~ say "here\n";
#~ $sequence = uc $_;
#~ #$sequence = ~ s/-/N/gi;
#~ print $FH_TOUT $_;
#~ }
#~ else {
#~ print $FH_TOUT $_;
#~ }
#~ }
#~ }
#~ print $FH_TOUT "\n";

#~ close $FH_IN;
#~ close $FH_TOUT;
#~ if ($flag eq "foobar"){
#~ my $my_command = "sort --output=$out_tbl_fn $out_tbl_fn";
#~ run($my_command);
#~ }
#~ #return $tblout;
#~ return 1;

#~ }

sub translate_2_protein ($genetic_code_fn, $cds_fn, $translated_tbl_fn)
{
    ### [<now>] Running translate_2_protein...
    ### $cds_fn
    ## 2016.12.12 my ($genetic_code_fn, $cds_fn, $translated_tbl_fn) = @_;

    #    print {*STDERR} "$genetic_code_fn in loop\n";

    my $frame = 0;

    #~ if ( !open(my $FH_FILEIN, "<", "$geneticcode" ) ) {
    #~ print "translate_2_protein: impossible to open genetic.code\n";
    #~ exit(1);
    #~ }

    my %gencodeh = ();
    open(my $FH_GENETIC_CODE, "<", "$genetic_code_fn")
      or croak "Can't open  $genetic_code_fn";
    while (<$FH_GENETIC_CODE>)
    {

        #        my $line = $_;

        #        my ($aa, $codon) = split(/\s+/, $line);
        my ($aa, $codon) = split;

        #print {*STDERR} "$codon\n";

        $gencodeh{$codon} = $aa;

        #print {*STDERR} "HERE: $gencodeh{$codon}\n";
    }
    close $FH_GENETIC_CODE;

    #if ( !open(my $FH_CDSIN, "<", "$cds" ) ) {
    #print "translate_2_protein: impossible to open $cds\n";
    #exit(1);
    #}

    open(my $FH_CDS_IN, "<", "$cds_fn")            or croak "err:  $cds_fn";
    open(my $FH_POUT,   ">", "$translated_tbl_fn") or croak "Failed here";

    #    print {*STDERR} "translating: $cds_fn \n";

    while (<$FH_CDS_IN>)
    {

        #        my $line = $_;

        #        my ($name, $seq) = split(/\s+/, $line);
        my ($name, $seq) = split;
        my $cds_seq_len = length($seq);

        print {$FH_POUT} "$name ";

        for (my $i = $frame ; $i < $cds_seq_len - 1 ; $i += 3)
        {
            my $triplet = substr($seq, $i, 3);
            if ($gencodeh{$triplet})
            {
                print {$FH_POUT} "$gencodeh{$triplet}";
            }
            else
            {
                print {$FH_POUT} "X";
            }

        }    #for
        print {$FH_POUT} "\n";
    }

    close $FH_CDS_IN;
    close $FH_POUT;
    ### translated_tbl_fn : $translated_tbl_fn
    ### [<now>] Finished translate_2_protein...
    return $translated_tbl_fn;

}

sub gff_2_geneidgff_mock ($my_input_nodots_gff, $species, $type)
{
    ## 2016.12.12 my ($my_input_nodots_gff, $species, $type) = @_;
    say "fake gff_2_geneidgff_mock\n";
    my $geneid_gff = $work_dir . $species . ${type} . ".geneid_gff";
    my $my_command = "cp $my_input_nodots_gff $geneid_gff";
    say "$my_command\n";
    run($my_command);

    my $geneid_sorted_gff_fn = $geneid_gff;
    return $geneid_sorted_gff_fn;
}
## CONVERT GFF2 TO GENEID GFF
# sub convert_GFF_2_geneidGFF {

#     my ( $my_input_nodots_gff, $species, $type ) = @_;
#     my %G;
#     my @G = ();

#     my $geneid_gff = $work_dir . $species . ${type} . ".geneid_gff";

#     open( my $FH_GFF,    "<", "$my_input_nodots_gff" ) or croak "Failed here";
#     open( my $FH_GFFOUT, ">", "$geneid_gff" )     or croak "Failed here";
#     while (<$FH_GFF>) {
#         my ( $c, @f, $id );
#         $c = ":";
#         $_ =~ s/\|//;
#         chomp;
#         @f  = split(/\s+/o);
#         $id = $f[8];           #seq name i.e. 7000000188934730
#         ( exists( $G{$id} ) ) || do {
#             $c = "#";
#             $G{$id} = [ @f[ 8, 0, 6 ], 0, @f[ 3, 4 ], [] ]
#               ;                # [7000000188934730 1.4 - 0 46549    46680 ]
#         };
#         push @{ $G{$id}[6] }, [ @f[ 3, 4 ] ];
#         $G{$id}[3]++;
#         $G{$id}[4] > $f[3] && ( $G{$id}[4] = $f[3] );
#         $G{$id}[5] < $f[4] && ( $G{$id}[5] = $f[4] );
#     }

#     @G =
#       sort { $a->[1] cmp $b->[1] || $a->[4] <=> $b->[4] || $a->[5] <=> $b->[5] }
#       map { $G{$_} } keys %G;

#     foreach my $GN (@G) {
#         my (
#             $id,     $ctg, $str, $nex, $go,   $ge, $coords,
#             @coords, $ce,  $CE,  $cur, $feat, $c
#         );
#         ( $id, $ctg, $str, $nex, $go, $ge, $coords ) = @{$GN};

#        # print {*STDERR} Data::Dumper->Dump([ \@$coords ], [ qw/ *coords / ])."\n";
#         @coords = map { $_->[0], $_->[1] }
#           sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
#           map { $_ } @{$coords};

#         # print {*STDERR} Data::Dumper->Dump([ \@coords ], [ qw/ *Coords / ])."\n";
#         #print "# $id $ctg $str $nex $go $ge\n";
#         $ce = 0;
#         $CE = $nex * 2;
#         $c  = 1;
#         while ( $ce < $CE ) {

#             # $cur = ($str eq '-' ? $CE - $ce - 2 : $ce);
#             if ( $nex == 1 ) {
#                 $feat = "Single";
#             }
#             elsif ( $c == 1 ) {
#                 $feat = $str eq '-' ? "Terminal" : "First";
#             }
#             elsif ( $c == $nex ) {
#                 $feat = $str eq '-' ? "First" : "Terminal";
#             }
#             else {
#                 $feat = "Internal";
#             }
#             print {$FH_GFFOUT} join( "\t",
#                 $ctg, "$species", $feat, $coords[$ce], $coords[ ( $ce + 1 ) ],
#                 ".", $str, ".", $id )
#               . "\n";
#             $ce += 2;
#             $c++;
#         }
#     }
#     close $FH_GFF;
#     close $FH_GFFOUT;

# #~ my $X_geneid_sorted_gff = "";
# #~ open( my $FH_LOCID, "-|", "sort -s -k8,9 -k4,5n $geneid_gff") or croak "Failed here";
# #~ while (<$FH_LOCID>) {
# #~ $X_geneid_sorted_gff .= $_;
# #~ }
# #~ close $FH_LOCID;

# ##
#     my $geneid_sorted_gff_fn =
#       $work_dir . $species . $type . ".geneid.gff_sorted";
#     my $my_command = "sort -s -k8,9 -k4,5n $geneid_gff > $geneid_sorted_gff_fn";
#     run($my_command);

#     #~ open( my $FH_FOUT, ">", "$geneid_sorted_gff_fn" ) or croak "Failed here";
#     #~ print {$FH_FOUT} "$geneid_sorted_gff_fn";
#     #~ close $FH_FOUT;

#     return $geneid_sorted_gff_fn;

#     #exit(0);

# }

sub sorteval
{
    ## descending numerical
    $b->[7] <=> $a->[7]    ## reverse numerical
      || $b->[10] <=> $a->[10]
      || $b->[4] <=> $a->[4]
      ## ascending numerical
      || $a->[11] <=> $b->[11]
      || $a->[12] <=> $b->[12]

}

#DK_subs
sub num_of_lines_in_file ($input_fn)
{
    ## 2016.12.11 my $input_fn = $_[0];
    my $my_num_lines;
    $my_num_lines = capture("cat $input_fn | wc -l ");

    #assert No such file or directory #
    chomp $my_num_lines;
    $my_num_lines = int($my_num_lines);
    return $my_num_lines;
}

sub check_external_progs()
{
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
          submatrix.py
          submatrix_order0.py
          get_k_matrix.py
          multiple_annot2one.awk
          logratio_kmatrix.py
          logratio_zero_order.py
          preparetrimatrixstart4parameter.awk)

    );

    for my $prog_name (@progs_2_check)
    {
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

sub create_data_dirs (@data_dirs)
{

    #my $dir_name = shift(@_);
    say "\ninside create_data_dirs function\n";
    ## check if not a bug XXX
    foreach my $dir_name (@data_dirs)
    {
        say "Creating $dir_name directory\n";
        if (-d $dir_name)
        {
            say "$dir_name exists. Purge and recreate\n";
            rmtree(["dir_name"]);
        }
        else
        {
            #say "$dir_name does not exist!\n";
            my $my_command = "mkdir -p $dir_name";
            run($my_command);
        }
    }
    return 1;
}

#sub write_sizes_from_tbl_fn {
#
#    my $input_tbl_fn = $_[0];
#    print "calc size for  $input_tbl_fn \n";
#
#    open( my $FH_input, '<', $input_tbl_fn )
#      or croak "Could not open file $input_tbl_fn' $OS_ERROR";
#
#    while ( my $line = <$FH_input> ) {
#        chomp $line;
#
#        #print $line
#        my @f = split / /, $line;
#        my $name = $f[0];
#
#        #print "$name\t"
#        my $seq       = $f[1];
#        my $lengfasta = length($seq);
#
#        #print "$name = $lengfasta\t"
#        my $tmp_out_fn = "$fastas_dir/$name" . "_len";
#
#        #print  " $tmp_out_fn \t"
#        open( my $FH_out, '>', $tmp_out_fn )
#          or croak "Could not open file '$tmp_out_fn' $OS_ERROR";
#        print {$FH_out} "$name $lengfasta\n";
#        close $FH_out;
#
#    }
#    close $FH_input;
#    return 1;
#}    #end write_sizes_from_tbl_fn

sub get_background_kmers ($kmer_len, $input_fas_fn, $contigs_all_tbl,
                          $num_seqs, $backgrnd_kmers_tbl)
{
    ## 206.12.13a {

    #~ my ($kmer_len, $input_fas_fn, $contigs_all_tbl, $num_seqs, $backgrnd_tbl) =
    #~ @_;
    my $backgrnd_debug_flag;
    $backgrnd_debug_flag = 1;
    if ($backgrnd_debug_flag)
    {
        ## TODO using python script to generate 1M 60mers to save time
        #just for non-random results
        print "\npremade background\n";
        my $my_command =
          "zcat ./test_data/pmar_1M_60mers.tbl.gz > $backgrnd_kmers_tbl";
        run($my_command);
        my $backgrnd_seq_freq_fn =
          "/tmp/workdir_00_gtrain/pmar01_background.freq";
        run("./bin/frequency.pypy 1 $backgrnd_kmers_tbl >  $backgrnd_seq_freq_fn"
        );
    }

    #    else  {
    #
    #        #my $input_fas_fn = $f;
    #        #my $contigs_all_tbl = "";
    #        my $countlines         = 0;
    #        my $contigs_all_number = num_of_lines_in_file($input_fas_fn);
    #
    #        #my $totalseqs  = `egrep -c \"^>\" $input_fas_fn`;
    #        #chomp $totalseqs;
    #
    #        open( my $FH_IN, "<", "$contigs_all_tbl" ) or croak "Failed here";
    #        my @tabular = ();
    #        while (<$FH_IN>) {
    #            push @tabular, $_;
    #        }
    #        close $FH_IN;
    #
    #        print {*STDERR}
    #"\nThe total number of genomic sequences (in $input_fas_fn) is $contigs_all_number\n";
    #
    #        if ( $contigs_all_number >= 1 ) {
    #
    #            #print {*STDERR} "in totalseqs if";
    #            my $seq    = "";
    #            my $sublen = 0;
    #            foreach my $line (@tabular) {
    #                chomp $line;
    #                my @f = split " ", $line;
    #                $seq .= $f[1];
    #                $sublen += length( $f[1] );
    #                $countlines++;
    #                if ( $sublen >= ( $num_seqs * $kmer_len + 1 ) ) {
    #                    last;
    #                }
    #            }
    #
    #            #print {*STDERR} "$seq";
    #            my $len = length($seq);
    #
    #            # chomp $totalseqs;
    #            print {*STDERR}
    #"The length of the new fasta file is $len\n(concatenation of $countlines fastas (out of $contigs_all_number))\n";
    #            open( my $FH_BACKGRND, ">", "$backgrnd_tbl" )
    #              or croak "Failed here";
    #            my $row = 1;
    #            print {*STDERR}
    #"\nObtain background ($num_seqs seqs) from fasta of $kmer_len nucleotides \n";
    #            for ( my $n = 1 ; $n <= $num_seqs ; $n++ ) {
    #                my $kmer = "N";
    #                while ( $kmer =~ m/[^ACGTacgt]/ ) {
    #                    my $rand = int( rand( $len - $kmer_len ) + 1 );
    #                    $kmer = substr( $seq, $rand, $kmer_len );
    #                }
    #
    #                #print {*STDERR}  $n."..";
    #                print {$FH_BACKGRND} $row++ . "\t$kmer\n";
    #            }
    #            print {*STDERR} "FINISHED OBTAINING BACKGROUND INFO\n";
    #            close $FH_BACKGRND;
    #        }
    #    }
    #

    return 1;

}    # END get_background_kmers function

# sub compute_sites_pictogram {
# #

#     my $my_command;
#     my $FH_FOUT;
# #
# ## get donor site statistics
# #
#     my $order = "0";
#     my $sites_number;
#     $sites_number = num_of_lines_in_file($out_donor_tbl);

#     my $don_offset =
#       $bases_offset;   #position before intron (last of exon (31) -1 for offset)

#     if ( $sites_number > $train_sites_cutoff ) {
#         $order = "1";
#     }
#     elsif ( $sites_number <= $train_sites_cutoff ) {
#         $order = "0";
#     }

#     print {*STDERR}
# "\nThere are $sites_number donor sites, enough for a matrix of order $order, prior offset: $don_offset $out_donor_tbl \n";

#     my (
#         $donor_matrix, $prof_len_don, $fxd_don_offset,
#         $donor_start,  $donor_end
#       )
#       = get_K_matrix( $out_donor_tbl, $backgrnd_kmers_fn, $order, $don_offset,
#         1, 0, 0, 0, 0, 0, 0 );
#     if (
#         !defined @{ $param->isocores }[0]->set_profile(
#             'Donor_profile', $prof_len_don, $fxd_don_offset, $pwm_cutoff,
#             $order, 0, 1, 0, 0, 0, 0, $donor_matrix
#         )
#       )
#     {
#         croak "error in setting profile\n";
#     }

# #print {*STDERR} "gawk '{print  substr(\$2,($donor_start-3),($prof_len_don+6))}' $out_donor_tbl\n";

#     $my_command =
# "gawk '{print  substr(\$2,($donor_start-3),($prof_len_don+6))}' $out_donor_tbl ";
#     my $donsub = capture($my_command);

#     my $donors_subprofile = $work_dir . $species . ".don.sub.profile";
#     open( $FH_FOUT, ">", "$donors_subprofile" ) or croak "Failed here";
#     print {$FH_FOUT} "$donsub";
#     close $FH_FOUT;

# # BUG pictogram
# #~ $my_command = "./bin/pictogram $donors_subprofile $plots_dir/donor_profile.pictogram -bits -land";
# #~ print "\n$my_command\n";
# #~ run($my_command);

# #
# # get acceptor site statistics
# #

#     $order = "0";

#     $sites_number = num_of_lines_in_file($out_acceptor_tbl);

#     #$sites_number = `wc -l $out_acceptor_tbl | gawk '{print \$1}'`;
#     #chomp $sites_number;
#     #$sites_number = int($sites_number);
#     print "sites_number in $out_acceptor_tbl: $sites_number\n";
#     my $acc_offset =
#       $bases_offset;   #position after intron (first of exon (31) -1 for offset)

#     if ( $sites_number > $train_sites_cutoff ) {

#         $order = "1";

#     }
#     elsif ( $sites_number <= $train_sites_cutoff ) {

#         $order = "0";
#     }

#     print {*STDERR}
# "\nThere are $sites_number acceptor sites, enough for a matrix of order $order, offset: $acc_offset \n";

#     my (
#         $acceptor_matrix, $prof_len_acc, $fxd_acc_offset,
#         $acceptor_start,  $acceptor_end
#       )
#       = get_K_matrix( $out_acceptor_tbl, $backgrnd_kmers_fn, $order,
#         $acc_offset, 0, 1, 0, 0, 0, 0, 0 );
#     if (
#         !defined @{ $param->isocores }[0]->set_profile(
#             'Acceptor_profile', $prof_len_acc, $fxd_acc_offset, $pwm_cutoff,
#             $order, 0, 1, 0, 0, 0, 0, $acceptor_matrix
#         )
#       )
#     {
#         croak "error in setting profile\n";
#     }

#     my $acceptors_subprofile = "";

#     $my_command =
# "gawk '{print  substr(\$2,($acceptor_start-3),($prof_len_acc+6))}' $out_acceptor_tbl ";
#     my $accsub = capture($my_command);

#     $acceptors_subprofile = $work_dir . $species . ".acc.sub.profile";
#     open( $FH_FOUT, ">", "$acceptors_subprofile" ) or croak "Failed here";
#     print {$FH_FOUT} "$acceptors_subprofile";
#     close $FH_FOUT;

# # BUG pictogram
# #~ $my_command = "./bin/pictogram $acceptors_subprofile $plots_dir/acceptor_profile.pictogram -bits -land";
# #~ print "\n$my_command\n";
# #~ run($my_command);

# #
# # get start site statistics
# #

#     $order       = "0";
#     $sites_number = num_of_lines_in_file($out_ATGx_tbl);

#     my $ATGx_offset =
#       $bases_offset;  #before first position of the exon (31)minus 1 for offset)

#     if ( $sites_number > $train_sites_markov_cutoff ) {

#         $order = "2";

#     }
#     elsif ( $sites_number <= $train_sites_markov_cutoff ) {

#         $order = "0";
#     }

#     print {*STDERR}
# "\nThere are $sites_number ATGx start sites, enough for a matrix of order $order, offset: $ATGx_offset \n";

#     my ( $start_matrix, $prof_len_sta, $fxd_ATGx_offset, $ATGx_start,
#         $ATGx_end ) =
#       get_K_matrix( $out_ATGx_tbl, $backgrnd_kmers_fn, $order, $ATGx_offset,
#         0, 0, 1, 0, 0, 0, 0 );

# # write to parameter file
#     if (
#         !defined @{ $param->isocores }[0]->set_profile(
#             'Start_profile', $prof_len_sta, $fxd_ATGx_offset, $pwm_cutoff,
#             $order, 0, 1, 0, 0, 0, 0, $start_matrix
#         )
#       )
#     {
#         croak "error in setting profile\n";
#     }
# #

#     my $stasub = "";

#     $my_command =
# "gawk '{print  substr(\$2,($ATGx_start-3),($prof_len_sta+6))}' $out_ATGx_tbl ";
#     $stasub = capture($my_command);

#     my $ATGx_subprofile = $work_dir . $species . ".ATGx.sub.profile";
#     open( $FH_FOUT, ">", "$ATGx_subprofile" ) or croak "Failed here";
#     print {$FH_FOUT} "$stasub";
#     close $FH_FOUT;

# ## BUG pictogram
# #~ $my_command = "./bin/pictogram $ATGx_subprofile $plots_dir/ATGx_profile.pictogram -bits -land";
# #~ print "\n$my_command\n";
# #~ run($my_command);

# ## end get start site statistics
#     return $donor_start, $donor_end, $acceptor_start, $acceptor_end,
#       $ATGx_start,
#       $ATGx_end;

# #
# }

sub compute_matrices_4sites ($my_input_table, $my_site_type)
{
    ## 2016.12.13a{
    ## hacked from compute_sites_pictogram
    ## no pictogram, just values

    ## my ($my_input_table, $my_site_type) = @_;

    my %types_hash = (
                      acceptor => 'Acceptor_profile',
                      donor    => 'Donor_profile',
                      ATGx     => 'Start_profile',
                      tr_stop  => 'Stop_profile',
                     );

    my $my_order     = "0";
    my $my_offset    = $bases_offset;
    my $sites_number = num_of_lines_in_file($my_input_table);

    if (($my_site_type eq 'acceptor') || ($my_site_type eq 'donor'))
    {
        if ($sites_number > $train_sites_cutoff)
        {
            $my_order = "1";    ### [<now>] $my_order Acc/Don site...
        }
    }
    elsif ($my_site_type eq 'ATGx')
    {
        if ($sites_number > $train_sites_markov_cutoff)
        {
            $my_order = "2";    ### [<now>] $my_order ATGx site...
        }
    }
    else
    {
        say "ERROR, unknown site type ";
    }

    my (
        $my_site_matrix, $my_site_profile_len, $my_site_offset,
        $my_site_start,  $my_site_end
       )
      = get_K_matrix($my_input_table, $backgrnd_kmers_fn, $my_order,
                     $my_offset, $my_site_type);
    ## TODO / BUG: possibly 0, 1, 0, 0, 0, 0, denotes a specific type of profile
    ## Nope, except branch profile all is fine
    ## the Stop profile is missing in the original sources
    #    ('Branch_point_profile', $branch_profile_len, $fxdbraoffset, -50, $order,     0, 1, 40, 10, 0, 0,  $branchmatrix)
    #    ('Start_profile',        $prof_len_sta, $fxdstaoffset, $cutoff, $order, 0, 1, 0, 0, 0, 0, $startmatrix)
    #    ('Acceptor_profile',     $prof_len_acc, $fxdaccoffset, $cutoff, $order, 0, 1, 0, 0, 0, 0, $acceptormatrix)
    #    ('Donor_profile',       $prof_len_don, $fxddonoffset, $cutoff, $order,  0, 1, 0, 0, 0, 0, $donormatrix)
    ## check inside the *.pm modules!
    if (
        !defined @{$param->isocores}[0]->set_profile(
                             $types_hash{$my_site_type},
                             $my_site_profile_len, $my_site_offset, $pwm_cutoff,
                             $my_order, 0, 1, 0, 0, 0, 0, $my_site_matrix
        )
       )
    {
        croak "error in setting profile of type: $my_site_type}\n";
    }

    return $my_site_start, $my_site_end;

}
