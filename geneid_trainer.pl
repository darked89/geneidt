#!/usr/bin/perl

## if run under perlbrew, use i.e.:
#!/$HOME/perl5/perlbrew/perls/perl-5.10.1/bin/perl

## checks & debugs modules
use Modern::Perl;
use English '-no_match_vars';

use feature "signatures";
no warnings "experimental::signatures";
use feature 'say';

use strict;
use warnings;
use autodie;
use diagnostics -verbose;
use sigtrap qw(stack-trace old-interface-signals);

use Carp qw(carp cluck croak confess);
use Carp::Always;
use Carp::Assert qw(assert);

use Smart::Comments;

#use Data::Dumper::Perltidy;
use Data::Dumper;

## common modules, used in part at the moment
#use Env qw(PATH, PERL5LIB);
use Cwd;

#use Getopt::Long;

use File::Path;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;

use YAML::XS qw(LoadFile);
use IPC::System::Simple qw(run system capture EXIT_ANY);
use Readonly;

## geneid_trained modules
use lib '.';
use Geneid::Param;
use Geneid::Isocore;
use Geneid::geneid;
use Geneid::geneidCEGMA;

#~ use Geneid::foo;

#~ #Geneid::GeneModel->new());
#~ Geneid::foo::baba();
#~ say "OK";
#~ die("test");

#~ use Inline Python => <<'END_OF_PYTHON_CODE2';

#~ def subtract(x,y):
#~ return x - y

#~ END_OF_PYTHON_CODE2

#~ die("XXX");
## MAIN VARIABLES
## my $PROGRAM_NAME    = "geneid_trainer";
## my $PROGRAM_VERSION = "2016.12.15a";

my $PROGRAM_HOME = getcwd;
my $exec_path    = "$PROGRAM_HOME/bin/";

local $ENV;
$ENV{'PATH'} = $exec_path . ":" . $ENV{'PATH'};

our $conf = LoadFile('gtrain_perl_conf.yaml');
print Dumper($conf);

#~ say $conf->{train_cds_fas};

## die("checking YAML");

my $genetic_code = "./etc/genetic.code";

## no need to run anything if this fails
##check_external_progs();

## Move parts necessary for getting comand line args here
my $input_hash   = $conf->{input};
my $input_gff_fn = $input_hash->{input_gff};
my $input_fas_fn = $input_hash->{genome_fasta};
my $species      = $input_hash->{species};

#~ my $input_gff_fn = $conf->{input_gff};
#~ my $input_fas_fn = $conf->{genome_fasta};
#~ my $species      = $conf->{species};

### Got:
### $species
### $input_gff_fn
### $input_fas_fn
###

## die("checking YAML");
### [<now>] Starting script at <file>[<line>]...

## TOFIX: move outside to YAML?
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

Readonly::Hash my %profile_params => (
    ## EXON WEIGHT PARAMETER
    #~ ExWeightParam_ini   => -5.0,
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

my $work_dir = $conf->{input}->{work_dir};
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

my $genome_all_contigs_tbl = "";

my @evaluation = ();

my $train_set_gff               = "";
my $eval_set_gff                = "";
my $contigs_all_transcr_2cols   = "";
my $train_contigs_transcr_2cols = "";
my $eval_contigs_transcr_2cols  = "";
my $train_transcr_used_int      = 0;
my $backgrnd_freq_fn            = $conf->{backgrn_freq};
my $backgrnd_kmers_tbl          = $conf->{backgrn_tbl};

## Golden Path stuff
my $gp_evalcontig_fa      = "";
my $gp_evalcontig_gff     = "";
my $gp_evalcontig_len_int = 0;
my $gp_evalcontig_tbl     = "";

my $gp_traincontig_fa      = "";
my $gp_traincontig_gff     = "";
my $gp_traincontig_len_int = 0;
my $gp_traincontig_tbl     = "";

## Weights for profiles??
my $best_ExWeightParam = 0;
my $best_OlWeight      = 0;
my $best_Acc           = 0;
my $best_Min           = 0;

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

my $train_2cols_seq_locusid_fn = "";
my $eval_2cols_seq_locusid_fn  = "";
my $train_transcr_ids_fn       = "";

my $train_noncanon_accept_int = 0;
my $train_noncanon_donors_int = 0;
my $train_noncanon_ATGx_int   = 0;

my $transcr_all_number = 0;
my $train_transcr_num  = 0;
my $eval_transcr_num   = 0;
## CREATE BLANK PARAMETER FILE##
### [<now>] template parameter file at <file>[<line>]
my $param        = Geneid::Param->new($species);
my $new_param_fn = "$work_dir/$species.geneid.param";

## END CREATING A PARAMETER FILE REGARDLESS OF WHETHER THE TRAINING IS COMPLETE OR REDUCED

normal_run();

sub normal_run
{
    our $conf;
    my $FH_FOUT;
    my $my_command = "";

    init_param_file();

    ### get data from upstream Python/YAML ###
    ## preserve input GFF file...
    $my_command = "cp $input_gff_fn $input_nodots_gff";
    run($my_command);

    ## number of genomic sequences TOTAL
    my $total_genomic = $conf->{all_contigs_int};

    ## number of gene models TOTAL
    $transcr_all_number = $conf->{all_transc_int};

    ## transcript ids all fn
    my $transcr_all_ids_fn = $work_dir . $species . "_all_seqs.lst";
    $my_command = "cp '$conf->{all_transc}' $transcr_all_ids_fn";

    ### Obtain 2cols for all transcripts at <file>[<line>]...
    $contigs_all_transcr_2cols = $work_dir . $species . "_locus_id";
    $my_command = "cp '$conf->{all_cont_transc}' $contigs_all_transcr_2cols";
    run($my_command);

    ## TOFIX -> move to outside Python
    ##assert $transcr_all_number >= $train_loci_cutoff)

    ### TRAIN/EVAL sets from upstream
    $train_transcr_num = $conf->{train_transc_int};
    $eval_transcr_num  = $conf->{eval_transc_int};

    $train_transcr_used_int = $train_transcr_num;   ## no filtration inside Perl

    ### TRAIN
    ## 2cols
    $train_contigs_transcr_2cols =
      $work_dir . $species . "_train_setaside80.2cols";
    $my_command =
      "cp '$conf->{train_cont_transc}' $train_contigs_transcr_2cols";
    run($my_command);

    ### gff
    my $gff_4_training = "";
    $train_set_gff = $work_dir . $species . "_train_setaside80.gff";
    $my_command    = "cp '$conf->{train_gff}' $train_set_gff";
    run($my_command);

    ### IDS
    $train_transcr_ids_fn = $work_dir . $species . "train_setaside80.ids";
    $my_command           = "cp '$conf->{train_transc}' $train_transcr_ids_fn";
    run($my_command);

    ### EVAL inputs
    ## 2cols
    $eval_contigs_transcr_2cols =
      $work_dir . $species . "_evaluation_setaside20.2cols";
    $my_command = "cp '$conf->{eval_cont_transc}' $eval_contigs_transcr_2cols";
    run($my_command);

    ## gff
    $eval_set_gff = $work_dir . $species . "_evaluation_setaside20.gff";
    $my_command   = "cp '$conf->{eval_gff}' $eval_set_gff";
    run($my_command);

    ### Not filtering premade inputs
    $train_cds_filtered_tbl       = $conf->{train_cds_tbl};
    $train_introns_filtered_tbl   = $conf->{train_intron_tbl};
    $train_inframestop_int        = 0;
    $train_locusid_filtered_2cols = $train_2cols_seq_locusid_fn;
    $train_filtered_gff           = $train_set_gff;

    #~ ### compability OLD, FIXME
    #~ $train_2cols_seq_locusid_fn = $train_set_gff;
    #~ gff_2_geneidgff_mock($train_set_gff, $species, ".train");
    #~ $eval_2cols_seq_locusid_fn = $eval_set_gff;
    #~ gff_2_geneidgff_mock($eval_set_gff, $species, ".eval");

    ### OLD FAS/TBL
    $genome_all_contigs_tbl = $work_dir . $species . ".genomic_all_contigs.tbl";
    $my_command =
      "./bin/fas_to_tbl.py $input_fas_fn $genome_all_contigs_tbl $fastas_dir";
    run($my_command);
    ## BUG?
    run("sort -o $genome_all_contigs_tbl $genome_all_contigs_tbl");

    ### [<now>] Getting premade sites tables/freq  <file>[<line>]...

    ### [<now>] Get background kmers freq premade  <file>[<line>]...
    #~ my $backgrnd_freq_fn = $conf->{backgrn_freq}; ## ERRORS...

    ### [<now>] Donor/Acc/ATG tables at <file>[<line>]...

    $train_donor_tbl           = $conf->{train_donor_tbl};
    $train_noncanon_donors_int = 0;
    $train_acceptor_tbl        = $conf->{train_acc_tbl};
    $train_noncanon_accept_int = 0;
    $train_ATGx_tbl            = $conf->{train_ATG_tbl};
    $train_noncanon_ATGx_int   = 0;

    my $train_donor_num    = $conf->{train_donor_num};
    my $train_acceptor_num = $conf->{train_acc_num};
    my $train_ATGx_num     = $conf->{train_ATG_num};

    my ($donor_start, $donor_end) =
      compute_matrices_4sites($train_donor_tbl, $train_donor_num, 'donor');
    my ($acceptor_start, $acceptor_end) =
      compute_matrices_4sites($train_acceptor_tbl, $train_acceptor_num,
                              'acceptor');
    my ($ATGx_start, $ATGx_end) =
      compute_matrices_4sites($train_ATGx_tbl, $train_ATGx_num, 'ATGx');

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

    ### PARAM settings
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

    my $train_introns_filtered_tbl = $conf->{train_intron_tbl};

    ##my @introns_len_list = seq_sizes_calc($train_introns_filtered_tbl);
    my ($intr_mean, $intr_stdev, $intr_short, $intr_long) =
      mean_stdev($train_introns_filtered_tbl);

    ### "INTRONS mean, ST

    ## $intron_short_int = 21;
    ## $intron_long_int  = 6384;
    ## => getting worse parameter

    my ($intron_short_int, $intron_long_int) =
      minmax_intron_special($intr_mean, $intr_stdev, $intr_short, $intr_long);
    my $intergenic_min = 200;
    my $intergenic_max = 'Infinity';

    ### PARAM INTRONS / INTERGEN
    ## MOVE UP?
    ## use shortest and longest intron lengths in gene model of parameter file
    $param->geneModel->intronRange($intron_short_int, $intron_long_int);
    $param->geneModel->intergenicRange($intergenic_min, $intergenic_max);
    ####################################################################

    ### [<now>] Write preliminary non-optimized parameter file at <file>[<line>]...

    $new_param_fn = "$work_dir/$species.geneid.param";
    $param->writeParam($new_param_fn);

    ### [<now>] Optimizing new parameter file at <file>[<line>]...

    say "new gff2gp location DEBUG";
    ### [<now>] Convert gff to gp (golden-path-like)format <file>[<line>]...
    ### artificial contig - concatenated sequences - approx. 800bp linkers...

    (
     $gp_traincontig_gff, $gp_traincontig_fa,
     $gp_traincontig_tbl, $gp_traincontig_len_int
    )
      = @{process_seqs_4opty($train_filtered_gff, ".train", 1)};

    ### [<now>] at <file>[<line>]...

    $eval_filtered_gff = $conf->{eval_gff};
    (
     $gp_evalcontig_gff, $gp_evalcontig_fa,
     $gp_evalcontig_tbl, $gp_evalcontig_len_int
    )
      = @{process_seqs_4opty($eval_filtered_gff, ".eval", 1)};

    ## OPTIMIZATION FUNCTION NO BRANCH
    my $array_ref         = "";
    my $OligoWeight_ini   = $profile_params{'OligoWeight_ini'};
    my $ExWeightParam_ini = $profile_params{'ExWeightParam_ini'};

    ### [<now>] contig optimization at <file>[<line>]...

    @evaluation = @{
        parameter_optimize(
                           $gp_traincontig_fa,
                           $gp_traincontig_gff,
                           $new_param_fn,
                           %profile_params
                          )
                   };

    ($best_ExWeightParam, $best_OlWeight, $best_Acc, $best_Min, $array_ref) =
      @{get_opt_paramfile(\@evaluation)};

    my @evaluation_init = @{$array_ref};
    my @evaluation_test = ();

    ##
    ## EVALUATE PERFORMANCE OF NEW PARAMETER FILE ON TEST SET (IF PRESENT)
    ##

    my $param_opt_fn = "$species.geneid.optimized.param";

    open(my $FH_SOUT, ">", "$work_dir/$species.use_NOT_allseqs.log");

    @evaluation_test = @{
        parameter_evaluate(
                           $gp_evalcontig_fa, $gp_evalcontig_gff, $param_opt_fn,
                           $OligoWeight_ini,  $ExWeightParam_ini
                          )
                        };

    print {*STDERR}
      "\nPerformance of new optimized parameter file on test set:\n\n"
      . join("\t", @evaluation_init[2 .. $#evaluation_init]), "\n";

    print {*STDERR} join("\t", @evaluation_test), "\n\n";

    print {$FH_SOUT} join("\t", @evaluation_test), "\n\n";
    close $FH_SOUT;

    return 1;
} ## end normal run

#
## END OF MAIN PORTION OF SCRIPT
#

sub init_param_file()
{
    ###########################################
    ### PARAM INIT
    #set isochores to 1
    $param->numIsocores(1);
    $param->isocores([Geneid::Isocore->new()]);

    $param->geneModel(Geneid::GeneModel->new());
    $param->geneModel->useDefault;
    ###########################################
}

## FUNCTION TO OBTAIN MARKOV MODELS CORRESPONDING TO THE CODING POTENTIAL
sub derive_coding_potential ($in_cds_tbl_fn, $in_intron_tbl_fn)
{
    ### [<now>] Running derive_coding_potential at <file>[<line>]...

    my $markov_mod_A         = "";
    my $markov_mod_B         = "";
    my $total_codingbases    = $conf->{train_cds_len};
    my $total_noncodingbases = $conf->{train_introns_len};

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
      ## BUG: independent of condition markov_mod_A&B always = 5 & 4
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
    ### [<now>] Finished derive_coding_potential at <file>[<line>]...
    return [
            \@profile_init,        \@profile_trans, $total_codingbases,
            $total_noncodingbases, $markov_mod_A
           ];

}    #derive coding potential

## PROCESS SEQUENCES FUNCTION ( FLANKED GENE MODELS OBTAINED FOR OPTIMIZATION)
sub process_seqs_4opty ($in_gff, $opt_type, $contig_opt_flag)
{

    ### [<now>] Running process_seqs_4opty at <file>[<line>]...

    my $pso_gp_tbl  = "";
    my $gp_from_gff = "";
    my $pso_gp_gff  = "";

    ##my $my_command  = "./bin/gff2gp.pypy $in_gff | sort -k 1";
    my $my_command = "./bin/gff2gp.pypy $in_gff";
    $gp_from_gff = capture($my_command);

    my $pso_tmp_gp_from_gff =
      $work_dir . $species . $opt_type . $contig_opt_flag . ".gp";

    open(my $FH_FOUT, ">", "$pso_tmp_gp_from_gff");
    print {$FH_FOUT} "$gp_from_gff";
    close $FH_FOUT;

    ### [<now>] Before get_genes at <file>[<line>]...
    ### $fastas_dir, $pso_tmp_gp_from_gff, $work_dir , $contig_opt_flag \n";

    my $gp_Xgetgenes_tmp_pre_tbl =
      get_genes($fastas_dir, $pso_tmp_gp_from_gff, $work_dir);

    ### [<now>] After get_genes at <file>[<line>]...
    ### PRETBL AFTER GETGENES: $gp_Xgetgenes_tmp_pre_tbl
    ### Get sequences of 400-nt flanked sequences in tabular and gff formats\n";

    $my_command =
      "gawk 'BEGIN{b=\"x\"}{if (\$1!=b){printf \"\\n\%s \",\$1}printf \"\%s\",\$6;b=\$1}END{printf \"\\n\"}' $gp_Xgetgenes_tmp_pre_tbl | sort | sed '/^\$/d' - | gawk '{print \$1, toupper(\$2)}' - ";
    $pso_gp_tbl = capture($my_command);

    my $temp_gp_tbl =
      $work_dir . $species . $opt_type . $contig_opt_flag . ".gp.tbl";
    open($FH_FOUT, ">", "$temp_gp_tbl") or croak "Failed here";
    print {$FH_FOUT} "$pso_gp_tbl";
    close $FH_FOUT;

    $my_command =
      "gawk 'BEGIN{OFS=\"\\t\";pos=1;b=\"x\"}{if (\$1!=b){pos=1}; print \$1,\"annotations\",\$3,pos,pos+\$5-1,\"\.\",\"+\",\"\.\",\$1\$2; pos+=\$5;b=\$1 }' $gp_Xgetgenes_tmp_pre_tbl | egrep -v '(Intron|Utr)' - ";
    $pso_gp_gff = capture($my_command);

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

    my $opt_type_cds_contigs_total_bp =
        $work_dir
      . $species
      . $opt_type
      . $contig_opt_flag
      . ".gp_gene400flank_length";
    open($FH_FOUT, ">", "$opt_type_cds_contigs_total_bp")
      or croak "Failed here";
    print {$FH_FOUT} "$gp_seqs_lengths";
    close $FH_FOUT;

    $my_command =
      "./bin/gff2cds.awk source=\"annotations\" $temp_gp_gff | sort -k1,1 | join $opt_type_cds_contigs_total_bp - ";
    my $cds_gp = capture($my_command);

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

    my @gp_tabular = split(/\n/, $pso_gp_tbl);
    my $seq = "";
    foreach my $line (@gp_tabular)
    {
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

    $my_command =
      "./bin/multiple_annot2one.awk species=$species leng=$my_seq_len $temp_gp_cds_fn ";
    my $gp_contig_tmp = capture($my_command);

    my $temp_gff2gpcontig =
      $work_dir . $species . $opt_type . $contig_opt_flag . ".contig.gp.cds";

    open($FH_FOUT, ">", "$temp_gff2gpcontig") or croak "Failed here";
    print {$FH_FOUT} "$gp_contig_tmp";
    close $FH_FOUT;

    $my_command =
      "./bin/cds2gff.awk $temp_gff2gpcontig | gawk 'BEGIN{OFS=\"\\t\";}{if (NR==1){print \"$species\",\"annotations\",\"Sequence\",\"1\",$my_seq_len,\".\",\".\",\".\",\".\";print}else {print}}' - ";
    my $cds2gffcontig = capture($my_command);

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
    ### Finished process_seqs_4opty at <file>[<line>]...
    say "L1302";
    return [
            $temp_gp_cdsgff_contig_eval, $temp_fastagpcontig,
            $temp_tabulargpcontig,       $temp_seqlencontig
           ];

}    #processSequences optimization

## GETGENES FUNCTION: EXTRACT FLANKED SEQUENCES FROM GENE MODELS FOR LATER OPTIMIZATION
sub get_genes ($my_fastas_dir, $my_pso_tmp_gp_from_gff, $work_dir)
{
    ### [<now>] Running get_genes at <file>[<line>]...
    my $only_non_red = 0;
    my $gp_out_tbl   = "$work_dir/gp_exon_utr_intron.tbl";

    open(my $FH_REFGENE, "<", "$my_pso_tmp_gp_from_gff")
      or croak "Failed here";
    ### [<now>] Reading at <file>[<line>]...
    ###  $my_pso_tmp_gp_from_gff

    open(my $FH_OUT_tblgb, ">", "$gp_out_tbl") or croak "Failed here";
    while (<$FH_REFGENE>)
    {
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

            #~ die("after first f2t run...");

            #print {*STDERR} "FATOTBL: $my_fastas_dir"."$chro\n";
            open(my $FH_IN, "<", "$chrom_tmp_tbl") or croak "Failed here";

            #            my @tabular = ();
            my $sub_seq = "";
            while (<$FH_IN>)
            {
                my @columns_f = split;
                $sub_seq .= $columns_f[1];

            }
            close $FH_IN;

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

            my $genomic_len = length($genomic_tmp_seqX);
            my $cds_seq     = "";

            if ($genomic_len == 0)
            {
                print {*STDERR} "getgenes: gene of 0 length ($name), $call\n";
                next;
            }

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
    ### [<now>] Finished get_genes at <file>[<line>]...
    return $gp_out_tbl;

}    #getgenes

use Inline Python => <<'END_OF_PYTHON_CODE3';

def bit_score_graph(info_fn, info_thresh, offset):
    line_counter = 0
    accepted_positions = []
    if info_thresh  in [ 0.15, 0.04, 0.15]:
        #all ok
        for line in open(info_fn).readlines():
            line_counter += 1
            sl = line.split()
            matrix_position = int(sl[0])
            info_score      = float(sl[1]) 
            if info_score > info_thresh:
                accepted_positions.append(matrix_position)
        accepted_positions.sort()    
        start = accepted_positions[0]
        end   = accepted_positions[-1]
        if (start > offset - 1):
            start = offset - 1
        
        ## BUG??
        ## it should be one base after the GT/AG/ATG etc 2fix???
        if (end < offset + 1):
            end = offset + 1
        
    return start, end

END_OF_PYTHON_CODE3

## GETKMATRIX FUNCTION (Splice sites and Start ATG codon PWMs)
sub get_K_matrix ($true_kmers_tbl, $backgrnd_freq_fn, $order, $offset,
                  $matrix_type)
{

    ### [<now>] Running get_K_matrix at <file>[<line>]...
    ### [<now>] gKm1 $true_kmers_tbl, $backgrnd_freq_fn
    ### [<now>] gKm2  $order, $offset, $matrix_type
    ## my $original_offset = $offset;
    my @profile_array = ();
    ## $temp_infolog;
    my @orders  = (qw(hmm_0 hmm_2 hmm_3 hmm_4 hmm_5 hmm_6 hmm_7 hmm_8));
    my $ordname = $orders[$order];
    my $sort    = "sort -n";
    use Inline Python => <<'END_OF_BACKGR_MATRICES_CALC';
def background_calc(backgound_tab_fn):
	import os
	## assumption/temp hack
	matrix_len = 60 
	hmm_order = 1
	backgrnd_freq_matrix_fn = 'foobar' #FIXME
	
	exe = './bin/get_k_matrix.pypy'
	command_A  = '%s %s %s %s' % (exe, hmm_order, matrix_len, backgrnd_kmers_tbl)  
	command_B  = '| sort > ' 
	command_C  = '%s' % (backgrnd_freq_matrix_fn)
	#os.system(command)
	command = ""
	#os.system(command)
	command = ""
	#os.system(command)
END_OF_BACKGR_MATRICES_CALC
    my ($donor, $acceptor_mtype, $ATGx) = (0, 0, 0);

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
        $ATGx = 1;    ### [<now>] $matrix_type ATGx at <file>[<line>]...
    }
    else
    {
        croak "sub_get_K_matrix: wrong matrix type\n";
    }

    if ($order > 1)
    {
        $sort = "sort -k1,1n -k2,2";
    }

    #    my @info = ($offset-1,$offset+1);
    my $profile_len = 0;

    #~ ## BUG?
    my $true_seq_name = $true_kmers_tbl;
    $true_seq_name =~ s/\.tbl$//;
    ### true_seq_name $true_seq_name;

    #~ ## Open true sequences
    #~ #    print {*STDERR} "$true_kmers_tbl (true)\n";
    #~ open(my $FH_TRUE_SEQ, "<", "$true_kmers_tbl") or croak "Failed here";
    #~ $_ = <$FH_TRUE_SEQ>;
    #~ my @columns_t = split;
    #~ my $len       = length($columns_t[1]);
    #~ close $FH_TRUE_SEQ;

    ## mock: not checking the sequence lenght in matrices.
    ## this is done upstream by py script
    my $len  = 60;
    my $len2 = 60;

    #    die "$len != $len2\n" if $len != $len2;
    my $true_seq_freq_fn = $work_dir . basename($true_seq_name) . ".freq";

    #my $subtracted_true_false_freq_fn = $work_dir . basename($true_seq_name) . "_" . basename($backgrnd_seq_name).freq_subtr";
    my $my_freq_subtract_fn =
      $work_dir . basename($true_seq_name) . "_" . "backgrnd" . ".information";

    run("./bin/frequency.pypy 1 $true_kmers_tbl  >  $true_seq_freq_fn");

    ### [<now>] not running background freq again at <file>[<line>]...

    my $my_command_A =
      "./bin/information.pypy  $true_seq_freq_fn $backgrnd_freq_fn ";

    my $my_command_C = " > $my_freq_subtract_fn ";

    #my $my_command = $my_command_A . $my_command_B;
    my $my_command = $my_command_A . $my_command_C;

    say "\n $my_command \n";
    run($my_command);

    ## True_False req_matrix_fn
    my $my_true_freq_matrix_fn =
      $work_dir . basename($true_seq_name) . "_" . "$ordname.matrix";
    ## False_True req_matrix_fn
    my $my_backgrnd_freq_matrix_fn =
      $work_dir . "backgrnd" . "_" . "$ordname.matrix";

    ## True logratio req_matrix_fn
    my $my_T_generic_logratio_freq_matrix_fn =
      $work_dir . basename($true_seq_name) . "$matrix_type.log.$ordname.matrix";

    if (!$order)
    {
        ### HMM generic (no order?)
        $my_command =
          ##"gawk -f ./bin/logratio_zero_order.awk $backgrnd_freq_fn $true_seq_freq_fn > $my_T_generic_logratio_freq_matrix_fn";
          "./bin/logratio_zero_order.pypy $backgrnd_freq_fn $true_seq_freq_fn > $my_T_generic_logratio_freq_matrix_fn";
        say "\n $my_command \n";
        run($my_command);

    }
    elsif ($order == 1)
    {
        ### HMM $order 1
        $my_command =
          "./bin/get_k_matrix.pypy $order $len $true_kmers_tbl | $sort > $my_true_freq_matrix_fn";

        #~ "./py_code/get_k_matrix_2+kmers.py $order $len $true_kmers_tbl | $sort > $my_true_freq_matrix_fn";
        say "\n $my_command \n";
        run($my_command);

        $my_command =

          "./bin/get_k_matrix.pypy $order $len2 $backgrnd_kmers_tbl | $sort > $my_backgrnd_freq_matrix_fn ";

        #~ "./py_code/get_k_matrix_2+kmers.py $order $len2 $backgrnd_kmers_tbl | $sort > $my_backgrnd_freq_matrix_fn";

        say "\n $my_command \n";
        run($my_command);

        ##$my_backgrnd_freq_matrix_fn = $conf->{backgrnd_hmm2_matrix};

        $my_command =
          "./bin/logratio_kmatrix.py $my_backgrnd_freq_matrix_fn $my_true_freq_matrix_fn > $my_T_generic_logratio_freq_matrix_fn ";
        say "\n $my_command \n";
        run($my_command);

    }
    else
    {
        ### HMM_2+ $order >
        $my_command =
          "bin/get_k_matrix.pypy $order $len $true_kmers_tbl | $sort > $my_true_freq_matrix_fn";

        #~ "./py_code/get_k_matrix_2+kmers.py $order $len $true_kmers_tbl | $sort > $my_true_freq_matrix_fn";

        say "\n $my_command \n";
        run($my_command);

        $my_command =

          "./bin/get_k_matrix.pypy $order $len2 $backgrnd_kmers_tbl | $sort > $my_backgrnd_freq_matrix_fn ";

        #~ "./py_code/get_k_matrix_2+kmers.py $order $len2 $true_kmers_tbl | $sort > $my_true_freq_matrix_fn";

        say "\n $my_command \n";
        run($my_command);

        $my_command =
          "./bin/logratio_kmatrix.py $my_backgrnd_freq_matrix_fn $my_true_freq_matrix_fn > $my_T_generic_logratio_freq_matrix_fn ";
        say "\n $my_command \n";
        run($my_command);

    }
    ### [<now>]
    ### $matrix_type, $my_info_thresholds{$matrix_type}\n";

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
        ### [<now>] Running donor order_oneplus at <file>[<line>]...
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
        ### [<now>] Running acceptor order_one_plus at <file>[<line>]...
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
        ### [<now>] Running ATGx order_two at <file>[<line>]...
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
    {    ### [<now>] Running REMAINING CASES at <file>[<line>]...

        # print {*STDERR} "$path/submatrix_order0.awk $start $end $true_seq_name-log.$ordname-matrix\n";
        $my_command =
          "./bin/submatrix_order0.py $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_matrix_4param_fn";
        say $my_command;
        run($my_command);

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

}    ##END GETKMATRIX FUNCTION (Splice site an Start codon PWMs)

## sub numerically { $a <=> $b }

#Optimize parameter file
sub parameter_optimize ($gp_fasta, $gp_gff_fn, $new_param_fn, %profile_params)
{

    my $ExWeightParam_ini   = $profile_params{'ExWeightParam_ini'};
    my $ExWeightParam_step  = $profile_params{'ExWeightParam_step'};
    my $ExWeightParam_limit = $profile_params{'ExWeightParam_limit'};
    ## EXON/OLIGO FACTOR PARAMETER
    my $OligoWeight_ini   = $profile_params{'OligoWeight_ini'};
    my $OligoWeight_step  = $profile_params{'OligoWeight_step'};
    my $OligoWeight_limit = $profile_params{'OligoWeight_limit'};

    my @evaluation_total    = ();
    my $myExWeightParam_tmp = $ExWeightParam_ini;
    my $myOlWeight_tmp      = $OligoWeight_ini;

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

            run("cat $temp_evalout_B_fn");
            my @evaluation_output;
            open(my $FH_IN, "<", "$temp_evalout_B_fn") or croak "Failed here";
            while (<$FH_IN>)
            {
                #                @evaluation_output = split " ";
                @evaluation_output = split;

                #@evaluation_output = split " ", $_;
            }
            close $FH_IN;

            #~ print Dump @evaluation_output;
            #~ #

            #~ my @evaluation_output = split " ",

            #~ #` ./bin/evaluation -sta $temp_geneid_pred_gff_fn $gp_gff_fn | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $IoWF, $I_ExWeight_F, \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}' `;

            push(@evaluation_total, \@evaluation_output);

        }    #end for_#2
        $myOlWeight_tmp = $OligoWeight_ini;
        ###print {*STDERR} "\n";

    }    #end for_#1

    return \@evaluation_total;

}
## end sub optimize parameter file

sub get_opt_paramfile ( $eval_array )
{
    ### [<now>] Running get_opt_paramfile at <file>[<line>]...
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

    open(my $FH_SOUT, ">", "$work_dir/$species.get_opt_paramfile.log")
      or croak "Failed here";

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
    my $optimized_geneid_param_fn = "$species.geneid.optimized.param";
    $param->writeParam($optimized_geneid_param_fn);

    close $FH_SOUT;
    ### [<now>] Finished get_opt_paramfile at <file>[<line>]...
    ### $optimized_geneid_param_fn

    return [$best_ExWeightParam, $best_OlWeight, 0, 0, \@evaluation_init];

}

sub parameter_evaluate ($gp_fasta, $gp_gff_fn, $new_param_fn, $OligoWeight_ini,
                        $ExWeightParam_ini)
{
    ### [<now>] Running  parameter_evaluate at <file>[<line>]...

    my $my_command;
    my $geneid_test_predict_gff_fn =
      "$work_dir/geneid_test_predictions." . basename($new_param_fn) . ".gff";

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
    ### [<now>] Finished sub  parameter_evaluate at <file>[<line>]...
    return \@evaluation_test;

}    # evaluate parameter function

#~ sub translate_2_protein ($genetic_code_fn, $cds_fn, $translated_tbl_fn)
#~ {
#~ ### [<now>] Running translate_2_protein...
#~ ### $cds_fn
#~ ## 2016.12.12 my ($genetic_code_fn, $cds_fn, $translated_tbl_fn) = @_;

#~ my $frame    = 0;
#~ my %gencodeh = ();

#~ open(my $FH_GENETIC_CODE, "<", "$genetic_code_fn")
#~ or croak "Can't open  $genetic_code_fn";
#~ while (<$FH_GENETIC_CODE>)
#~ {
#~ my ($aa, $codon) = split;
#~ $gencodeh{$codon} = $aa;
#~ }
#~ close $FH_GENETIC_CODE;

#~ open(my $FH_CDS_IN, "<", "$cds_fn")            or croak "err:  $cds_fn";
#~ open(my $FH_POUT,   ">", "$translated_tbl_fn") or croak "Failed here";

#~ #    print {*STDERR} "translating: $cds_fn \n";

#~ while (<$FH_CDS_IN>)
#~ {
#~ my ($name, $seq) = split;
#~ my $cds_seq_len = length($seq);

#~ print {$FH_POUT} "$name ";

#~ for (my $i = $frame ; $i < $cds_seq_len - 1 ; $i += 3)
#~ {
#~ my $triplet = substr($seq, $i, 3);
#~ if ($gencodeh{$triplet})
#~ {
#~ print {$FH_POUT} "$gencodeh{$triplet}";
#~ }
#~ else
#~ {
#~ print {$FH_POUT} "X";
#~ }

#~ }    #for
#~ print {$FH_POUT} "\n";
#~ }

#~ close $FH_CDS_IN;
#~ close $FH_POUT;
#~ ### translated_tbl_fn : $translated_tbl_fn
#~ ### [<now>] Finished translate_2_protein...
#~ return $translated_tbl_fn;

#~ }

#~ sub gff_2_geneidgff_mock ($my_input_nodots_gff, $species, $type)
#~ {
#~ say "fake gff_2_geneidgff_mock\n";
#~ my $geneid_gff = $work_dir . $species . ${type} . ".geneid_gff";
#~ my $my_command = "cp $my_input_nodots_gff $geneid_gff";
#~ say "$my_command\n";
#~ run($my_command);

#~ my $geneid_sorted_gff_fn = $geneid_gff;
#~ return $geneid_sorted_gff_fn;
#~ }

#DK_subs

#~ sub check_external_progs()
#~ {
#~ ## Checking if the external programs are in the path.
#~ ## C, awk, python programs
#~ #my $prog_name;
#~ my $my_command;
#~ my @progs_2_check = (
#~ qw(bash
#~ gawk
#~ egrep
#~ sed
#~ geneid
#~ ssgff
#~ shuf
#~ pictogram
#~ gff2gp.awk
#~ cds2gff.awk
#~ frequency.py
#~ information.py
#~ submatrix.py
#~ submatrix_order0.py
#~ get_k_matrix.py
#~ multiple_annot2one.awk
#~ logratio_kmatrix.py
#~ logratio_zero_order.py
#~ preparetrimatrixstart4parameter.awk)

#~ );

#~ for my $prog_name (@progs_2_check)
#~ {
#~ $my_command = "which $prog_name > /dev/null";
#~ run($my_command);

#~ #or croak "\n $prog_name not available\n";
#~ }
#~ return 1;
#~ }

sub compute_matrices_4sites ($my_input_table, $sites_number, $my_site_type)
{
    ### [<now>] running compute_matrices_4sites
    ### sites_number: $sites_number
    #~ my $sites_number = num_of_lines_in_file($my_input_table);

    my $my_offset = $bases_offset;

    my $my_order =
      order_site_select($sites_number, $my_site_type, $train_sites_cutoff,
                        $train_sites_markov_cutoff);
    ### L2148 $my_input_table, $my_order, $my_offset, $my_site_type

    use Inline Python => <<'END_OF_COMPUTE_MATRICES_4SITES';

def order_site_select(sites_number, site_type, train_sites_cutoff, train_sites_markov_cutoff):
    ##sites_number = num_of_lines_in_file($my_input_table);
    #my_order = 0 #default
    encoding = 'utf-8'
    print("inside order_site_select")
    
    # fixing b-string:
    site_type = site_type.decode(encoding)
    if site_type in ('acceptor','donor'):
       if sites_number > train_sites_cutoff:
           my_order = 1
       else:
           my_order = 1
    elif site_type ==   'ATGx':
        if sites_number > train_sites_markov_cutoff:
            my_order = 2
        else:
            my_order = 0
    else:
        print( 'ERROR, unknown site',  site_type)      
    return my_order

END_OF_COMPUTE_MATRICES_4SITES

    my %types_hash = (
                      acceptor => 'Acceptor_profile',
                      donor    => 'Donor_profile',
                      ATGx     => 'Start_profile',
                      tr_stop  => 'Stop_profile',
                     );

    my (
        $my_site_matrix, $my_site_profile_len, $my_site_offset,
        $my_site_start,  $my_site_end
       )
      = get_K_matrix($my_input_table, $backgrnd_freq_fn, $my_order,
                     $my_offset, $my_site_type);
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

# back to simple subs
#

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

##  mean_stdev...
use Inline Python => <<'END_STATS';
from statistics import mean, stdev

def seq_sizes_calc(input_seq_table):
    sizes_list = []
    for line in open(input_seq_table).readlines():
        sl = line.split()
        seq_size = len(sl[1])
        sizes_list.append(seq_size)
        
    sizes_list.sort()
    print(sizes_list[:10])        
    return sizes_list
    
def mean_stdev(input_seq_table):
    input_list = seq_sizes_calc(input_seq_table)
    short_intr = input_list[0]
    long_intr  = input_list[-1]
    input_list = [float(x) for x in input_list]
    list_mean = mean(input_list)
    list_stdev = stdev(input_list)
    return list_mean, list_stdev, short_intr, long_intr
    
    
def minmax_intron_special( list_mean, list_stdev, short_intr, long_intr):    
    #~ slice1 = input_list[:5];
    #~ slice2 = input_list[-5:]
    
    short_intr = short_intr * 0.75
    if short_intr > 40:
        short_intr = 40 
    
    long_intr = list_mean + list_stdev * 3
    if long_intr > 100000:
        long_intr = 100000
    return short_intr, long_intr
        
END_STATS

sub create_data_dirs (@data_dirs)
{
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

sub tbl_2_fasta ($in_tbl_fn, $fa_out_fn)
{
    ### [<now>] Running tbl_2_fasta at <file>[<line>]...
    ### in: $in_tbl_fn
    ### out: $fa_out_fn

    open(my $FH_IN,   "<", "$in_tbl_fn") or croak "Failed here";
    open(my $FH_FOUT, ">", "$fa_out_fn") or croak "Failed here";
    while (<$FH_IN>)
    {
        chomp;

        #my ( $n, $s ) = split( /\s+/, $_ );
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
    ### [<now>] Finished tbl_2_fasta at <file>[<line>]...
    return $fa_out_fn;
}

sub tbl_2_single_fastas ($in_tbl_fn, $out_fas_dir)
{
### [<now>] Running tbl_2_single_fastas at <file>[<line>]...

    open(my $FH_IN_TBL, "<", "$in_tbl_fn") or croak "Failed here";

    #    print {*STDERR} "##tbl_2_single_fastas $in_tbl_fn\n";
    while (<$FH_IN_TBL>)
    {
        #        my $input_line = $_;
        #        chomp($input_line);

        #my ( $seq_name, $seq ) = split( /\s+/o, $_ );
        #        my ($seq_name, $seq) = split(/\s+/o, $input_line);
        my ($seq_name, $seq) = split;

        open(my $FH_FOUT_FASTA, ">", "${out_fas_dir}" . "$seq_name");

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
### [<now>] Finished tbl_2_single_fastas at <file>[<line>]...
    #close $FH_FOUT_FASTA;
    return 1;
}

sub fasta_2_tbl ($in_fa, $out_tbl_fn)
{
### [<now>] Running fasta_2_tbl at <file>[<line>]...
### $in_fa, $out_tbl_fn

    open(my $FH_IN,   "<", "$in_fa");
    open(my $FH_TOUT, ">", "$out_tbl_fn");

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
### [<now>] Finished fasta_2_tbl at <file>[<line>]...
    return $out_tbl_fn;
}

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

