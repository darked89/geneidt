#!/usr/bin/perl

## if run under perlbrew, use i.e.:
#!/$HOME/perl5/perlbrew/perls/perl-5.10.1/bin/perl

## checks & debugs modules
use Modern::Perl;
use English '-no_match_vars';

use feature "signatures";
no warnings "experimental::signatures";

use warnings::unused; # installs the check routine as 'once'

use warnings 'once';  # enables  the check routine

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
use Getopt::Long;

use File::Path;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;

use YAML::XS qw(LoadFile); #from dev_unstable

use IPC::System::Simple qw(run system capture EXIT_ANY);
use Readonly;

use YAML qw(Dump Bless);


use Function::Parameters qw(:strict);

## geneid_trained modules
use lib '.';
use Geneid::Param;
use Geneid::Isocore;
use Geneid::geneid;
use Geneid::geneidCEGMA;

my $PROGRAM_HOME = getcwd;

my $exec_path = "$PROGRAM_HOME/bin/";

local $ENV;
$ENV{'PATH'} = $exec_path . ":" . $ENV{'PATH'};

## UNUSED NOW my $genetic_code = "./etc/genetic.code";


## no need to run anything if this fails
check_external_progs();

## Move parts necessary for getting comand line args here
my $species      = "";
my $input_gff_fn = "";
my $input_fas_fn = "";
my $sout         = "-";
my $fix_fasta    = 0;

## Get arguments (command line)
GetOptions(
    'species:s'       => \$species,
    'gff:s'            => \$input_gff_fn,
    'fasta:s'         => \$input_fas_fn,
    'sout|statsout:s' => \$sout,
    'fix_fasta'        => \$fix_fasta

);
my $usage =
    "Usage: $PROGRAM_NAME -species H.sapiens -gff gff_name -fasta fasta_name -sout statsfile_out -fix_fasta 1";



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

### [<now>] Starting script at <loc>...


## end of getting ARGS
## set CONSTANTS
## Constant values. modify if needed
Readonly::Scalar my $pwm_cutoff     => -7;
Readonly::Scalar my $bases_offset   => 30;    #bases in front/after? a feature
## UNUSED HERE
## Readonly::Scalar my $train_fraction => 0.8;   #fraction of seq used for training
## Readonly::Scalar my $train_loci_cutoff         => 500;
Readonly::Scalar my $train_sites_cutoff         => 1_400;
## Readonly::Scalar my $train_sites_cutoff_alt  => 1_200; # changed in some part?
Readonly::Scalar my $train_sites_markov_cutoff  => 5_500;
Readonly::Scalar my $backgrnd_kmer_size        => 60;
Readonly::Scalar my $backgrnd_kmer_num         => 100_000;

#~ ( $totalcodingbases > 400000 && $totalnoncodingbases > 100000 )
Readonly::Scalar my $coding_bp_limit_A => 400_000;
Readonly::Scalar my $coding_bp_limit_B => 375_000;

#~ || ( $totalcodingbases > 375000 && $totalnoncodingbases > 150000 )
Readonly::Scalar my $non_coding_bp_limit_A => 100_000;
Readonly::Scalar my $non_coding_bp_limit_B => 150_000;
Readonly::Scalar my $non_coding_bp_limit_C => 35_000;

## UNUSED now Readonly::Scalar my $multi_total_noncodingbases => 25;
Readonly::Scalar my $coding_param_X => 0.25;
Readonly::Scalar my $intron_param_X => 10;

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
## UNUSED NOW
#~ Readonly::Hash my %canonical_const => (
        #~ Donor    => [31, 'GT'],
        #~ Acceptor => [28, 'AG'],
        #~ Start    => [30, 'ATG'],
        ## stop 2do
        ## branch 2do
##    );
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
my $work_dir     = "./workdir_gt_pre_20200319_001/";
my $new_param_fn = "$work_dir/$species.geneid.param";

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
    $cds_dir,   $geneid_dir, $results_dir
);

create_data_dirs(@data_dirs);

## PROGRAM SPECIFIC VARIABLES (unordered...)

my $input_nodots_gff = "$work_dir/input_gff_no_dots.gff";


my $genome_all_contigs_tbl = "";
my $backgrnd_kmers_fn      = "";

my @evaluation = ();

my $train_set_gff               = "";
my $eval_set_gff                = "";
my $contigs_all_transcr_2cols   = "";
my $train_contigs_transcr_2cols = "";
my $eval_contigs_transcr_2cols  = "";
my $train_transcr_used_int      = 0;

## Golden Path stuff
my $gp_evalcontig_fa      = "";
my $gp_evalcontig_gff      = "";
my $gp_evalcontig_len_int = 0;
my $gp_evalcontig_tbl     = "";


my $gp_traincontig_fa      = "";
my $gp_traincontig_gff      = "";
my $gp_traincontig_len_int = 0;
my $gp_traincontig_tbl     = "";


## Weights for profiles??
my $best_ExWeightParam = 0;
my $best_OlWeight      = 0;
my $best_Acc           = 0;
my $best_Min           = 0;


## my $X_geneid_sorted_gff = "";
my $seqs_eval_gff         = "";
my $train_inframestop_int = 0;
my $eval_inframestop_int  = 0;
## UNUSED 
## my $locus_id              = "";
##my $locus_id_new          = "";
my $intron_long_int       = 0;
my $intron_short_int      = 0;
my $intergenic_max        = 0;
my $intergenic_min        = 0;

## UNUSED now
## my $eval_cds_filtered_tbl       = "";
my $eval_filtered_gff           = "";
## UNUSED
##my $eval_introns_filtered_tbl   = "";
##my $eval_locusid_filtered_2cols = "";

my $train_acceptor_tbl           = "";
my $train_ATGx_tbl               = "";
my $train_cds_filtered_tbl        = "";
my $train_donor_tbl              = "";
my $train_filtered_gff           = "";
my $train_introns_filtered_tbl   = "";
my $train_locusid_filtered_2cols = "";

my $train_2cols_seq_locusid_fn = "";
my $eval_2cols_seq_locusid_fn  = "";
my $train_transcr_lst_fn       = "";

my $train_noncanon_accept_int  = 0;
my $train_noncanon_donors_int  = 0;
my $train_noncanon_ATGx_int    = 0;

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

#~ #-#  CREATE BLANK PARAMETER FILE...
#~ #-#  [<now>] template parameter file
#~ my $param        = Geneid::Param->new($species);
#~ my $new_param_fn = "$work_dir/$species.geneid.param";

#~ #set isochores to 1
#~ $param->numIsocores(1);
#~ $param->isocores([Geneid::Isocore->new()]);

#~ #-#  END CREATING A PARAMETER FILE
#~ #-#  [<now>] template parameter file at <file>[<line>]

my $param = init_param_file($species);

normal_run();

sub normal_run
{
    
    
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

    print {*STDERR} "\nObtain locus_id (list of genomic sequences / genes)\n";

    
    $contigs_all_transcr_2cols = "./precomputed/all_cont_transc.2col"; # @@@
        
    $transcr_all_number = 765; # @@@
    

    my $total_genomic = 765; # @@@
    
    #~ chomp $total_genomic;

    print {*STDERR}
        "\nThe gff file ($input_nodots_gff) contains a total of $total_genomic genomic sequences and $transcr_all_number gene models\n";

    #-#  get a list of transcripts TOTAL
    print {*STDERR} "\nObtain list of all transcripts\n\n";

    ## UNUSED
    ## my $transcr_all_list_fn = "./precomputed/pmar_all_transcr.ids";
    


    #if ($transcr_all_number >= $train_loci_cutoff)
    
    
    # FIXME 20200313 temp solution
    ## $transcr_all_number = 765;
    ## hack for if
    if ($transcr_all_number > 1)
    {

        ## $train_transcr_num = int($train_fraction * $transcr_all_number);
        $train_transcr_num = 612;
        print {*STDERR}
            "\nA subset of $train_transcr_num sequences (randomly chosen from the $transcr_all_number gene models) was used for training\n";
        ## DEBUG KEEP !!! shuf => random select
        ## head -$train_transcr_num just the first ones
        #my $my_command =           "shuf --head-count=$train_transcr_num $contigs_all_transcr_2cols | sort ";

        
        $train_contigs_transcr_2cols = "./precomputed/train_cont_transc.2col"; #@@@
        
        $train_set_gff = "precomputed/train.gff"; #@@@
        
        print {*STDERR}
            "\nThe new training gff file includes $train_transcr_used_int gene models (80% of total seqs)\n";

        print {*STDERR} "\nObtain list of training genes\n\n";
        
        $train_transcr_lst_fn = "./precomputed/pmar_train_transcr.ids"; # @@@
        
        $eval_contigs_transcr_2cols =   "./precomputed/eval_cont_transc.2col"; # @@@
        $seqs_eval_gff = 153; # @@@
        
        print {*STDERR}
            "The evaluation gff file includes $seqs_eval_gff gene models (20% of total seqs)\n\n";
        
        $eval_set_gff = "./precomputed/eval.gff"; # @@@
        

    }    # seqs > 500

    
    

        ## BUG => there is no need to convert gff to geneid format each time we run

        $train_2cols_seq_locusid_fn =
            gff_2_geneidgff_mock($train_set_gff, $species, ".train");
        print {*STDERR}
            "L575 : $train_2cols_seq_locusid_fn \t $train_contigs_transcr_2cols \n";

        ## FIXME 20200316: encode these in the yaml
        
            $train_cds_filtered_tbl        =  "./precomputed/pmar01.train.cds.tbl";
            $train_introns_filtered_tbl    =  "./precomputed/pmar01.train.introns.tbl";
            $train_locusid_filtered_2cols  =  "./precomputed/train_locusid_filtered.tbl";
            $train_filtered_gff             =  "./precomputed/train.gff";
            $train_inframestop_int        = 0; 
        

        $eval_2cols_seq_locusid_fn =
            gff_2_geneidgff_mock($eval_set_gff, $species, ".eval");

        print {*STDERR}
            "L533 tmp_locus_id_X_new:  $eval_2cols_seq_locusid_fn \t $eval_contigs_transcr_2cols\n";


            $eval_filtered_gff = "./precomputed/eval.gff"; #@@@
            $eval_inframestop_int = 0; # @@@
        

    ## extract and check splice sites and start codon. Use only canonical info #IN SEQUENCES USED IN TRAINING
    print {*STDERR}
        "L653 :  $train_filtered_gff \t $train_locusid_filtered_2cols  \n";
    
        $train_donor_tbl             = "./precomputed/pmar01.train.don_sites.tbl";  #@@@
        $train_noncanon_donors_int   = 0; #@@@
        $train_acceptor_tbl          = "./precomputed/pmar01.train.acc_sites.tbl";  #@@@
        $train_noncanon_accept_int   = 0; #@@@ 
        $train_ATGx_tbl              = "./precomputed/pmar01.train.start_sites.tbl";  #@@@ 
        $train_noncanon_ATGx_int     = 0; #@@@
     
    

 
    #-#  GET BACKGROUND SEQUENCES

    print
        "Obtaining $backgrnd_kmer_num background sequences of $backgrnd_kmer_size nucleotides each for estimating background frequencies of nucleotides\n";

    $backgrnd_kmers_fn = $work_dir . $species . "_background.tbl";
    get_background_kmers(
        $backgrnd_kmer_size,     
        $input_fas_fn,
        $genome_all_contigs_tbl, 
        $backgrnd_kmer_num,
        $backgrnd_kmers_fn
    );


    my ($donor_start, $donor_end) =
        compute_matrices_4sites($train_donor_tbl, 'donor');
    my ($acceptor_start, $acceptor_end) =
        compute_matrices_4sites($train_acceptor_tbl, 'acceptor');
    my ($ATGx_start, $ATGx_end) =
        compute_matrices_4sites($train_ATGx_tbl, 'ATGx');


    ### DERIVE INITIAL/TRANSITION MARKOV MODEL

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
    ## add markov matrices to the parameter file
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

            $train_introns_filtered_tbl,  $train_cds_filtered_tbl,
            $train_filtered_gff,           $train_inframestop_int,
            $eval_inframestop_int,       $train_transcr_used_int,
            $train_noncanon_donors_int,  $train_noncanon_accept_int,
            $train_noncanon_ATGx_int,    $markov_model,
            $total_coding,               $total_noncoding,
            $donor_start,                $donor_end,
            $acceptor_start,             $acceptor_end,
            $ATGx_start,                 $ATGx_end,
        );

    print {*STDERR}
        "\nshortest intron: $intron_short_int\nlongest intron: $intron_long_int\nminimum intergenic: $intergenic_min\nmaximum intergenic: $intergenic_max\n";

    ##
    #-#  WRITE PRELIMINARY NON-OPTIMIZED PARAMETER FILE
    #-#  "INTRONS mean, ST

    ## $intron_short_int = 21;
    ## $intron_long_int  = 6384;
    ## => getting worse parameter
    ### [<now>] Write preliminary non-optimized parameter file at <file>[<line>]...

    $new_param_fn = "$work_dir/$species.geneid.RAW.param";
    $param->writeParam($new_param_fn);


    ## # START prepare sequences for optimization BLOCK
    ## of newly developed parameter file (TRAIN)

    ### [<now>] TART prepare sequences for optimization at <file>[<line>]...

    print {*STDERR}
        "\nConvert gff to gp (golden-path-like)format (artificial contig - concatenated sequences - approx. 800 nt between sequences)\n";

    (
        $gp_traincontig_gff, $gp_traincontig_fa,
        $gp_traincontig_tbl, $gp_traincontig_len_int
    )
        = @{process_seqs_4opty($train_filtered_gff, ".train", 1)};
    print {*STDERR}
        "\nConvert gff to gp (golden-path-like)format (test set for evaluation of new parameter file -
        (artificial contig - concatenated sequences - approx. 800 nt between sequences)\n";
    (
        $gp_evalcontig_gff, $gp_evalcontig_fa,
        $gp_evalcontig_tbl, $gp_evalcontig_len_int
    )
        = @{process_seqs_4opty($eval_filtered_gff, ".eval", 1)};
    print {*STDERR} "DONE\n";
    
    ## # END prepare sequences for optimization BLOCK


    ### [<now>] Optimizing new parameter file at <file>[<line>]...

    ## OPTIMIZATION FUNCTION NO BRANCH
    my $array_ref         = "";
    my $OligoWeight_ini   = $profile_params{'OligoWeight_ini'};
    my $ExWeightParam_ini = $profile_params{'ExWeightParam_ini'};

    ### [<now>] contig optimization at <file>[<line>]...

    print {*STDERR} "\n DEBUG: CONTIG OPTIMISATION\n";

    @evaluation = @{
        parameter_optimize(
            $gp_traincontig_fa,
            $gp_traincontig_gff,
            $new_param_fn,
            %profile_params
        )
        };

    ($best_ExWeightParam, $best_OlWeight, $best_Acc, $best_Min, $array_ref) = @{
        BuildOptimizedParameterFile(\@evaluation)
        };

    my @evaluation_init = @{$array_ref};
    my @evaluation_test = ();

    ##
    #-#  EVALUATE PERFORMANCE OF NEW PARAMETER FILE ON TEST SET (IF PRESENT)
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
sub init_param_file($species)
{
### CREATE BLANK PARAMETER FILE...
### [<now>] template parameter file
my $param        = Geneid::Param->new($species);
## my $new_param_fn = "$work_dir/$species.geneid.param";

#set isochores to 1
$param->numIsocores(1);
$param->isocores([Geneid::Isocore->new()]);

### END CREATING A PARAMETER FILE
### [<now>] template parameter file at <file>[<line>]

    ### [<now>] running init_param_file at <file>[<line>]...

    ### PARAM INIT
    #set isochores to 1
    $param->numIsocores(1);
    $param->isocores([Geneid::Isocore->new()]);

    $param->geneModel(Geneid::GeneModel->new());
    $param->geneModel->useDefault;
    
    return $param;
}




sub derive_coding_potential ($train_cds_filtered_tbl, 
                             $train_introns_filtered_tbl)
{   #-#  FUNCTION TO OBTAIN MARKOV MODELS CORRESPONDING TO THE CODING POTENTIAL
    ### [<now>] Running derive_coding_potential at <file>[<line>]...
    ## my ($train_cds_filtered_tbl, $train_introns_filtered_tbl) = @_;

    my $markov_mod_A = "";
    my $markov_mod_B = "";

    my $my_command =
        "gawk '{ l=length(\$2); L+=l;} END{ print L;}' $train_cds_filtered_tbl";
    my $total_codingbases = capture($my_command);
    chomp $total_codingbases;

    $my_command =
        "gawk '{ l=length(\$2); L+=l;} END{ print L;}' $train_introns_filtered_tbl ";
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
            && $total_codingbases > (25 * $total_noncodingbases))
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

    open(my $FH_INTRONS_tbl, "<", "$train_introns_filtered_tbl") or croak "Failed here";
    my @intron_seqs = ();
    while (<$FH_INTRONS_tbl>)
    {
        my @columns_i = split;

        #print {*STDERR} "SECOND FIELD $i[2]";
        push @intron_seqs, $columns_i[1];
    }
    close $FH_INTRONS_tbl;

    open(my $FH_CDSes_tbl, "<", "$train_cds_filtered_tbl") or croak "Failed here";
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

}    #derive coding potential END

## PROCESS SEQUENCES FUNCTION ( FLANKED GENE MODELS OBTAINED FOR OPTIMIZATION)
sub process_seqs_4opty ($my_input_nodots_gff, $opt_type, $contig_opt_flag)
{
    ### [<now>] process_seqs_4opty at <file>[<line>]...
    ## my ($my_input_nodots_gff, $opt_type, $contig_opt_flag) = @_;

    #my $pso_out_tbl = "";
    my $pso_gp_tbl  = "";
    my $gp_from_gff = "";

    #~ my $gp_fasta    = ""; #unused
    my $pso_gp_gff = "";
    my $my_command = "";

    #my $work_dir;

    open(
        my $FH_LOCID, "-|",
        ##"./bin/gff2gp.py $my_input_nodots_gff | sort -k 1 " );
        "./bin/gff2gp.awk $my_input_nodots_gff | sort -k 1 "
    );

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

    ### [<now>] Before get_genes at <file>[<line>]...
    ### $fastas_dir, $pso_tmp_gp_from_gff, $work_dir , $contig_opt_flag \n";

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

    #  ` gawk '{print \$1,length(\$2)}' $temp_gp_tbl | sort -k1,1 `;    ##XX

    my $opt_type_cds_contigs_total_bp =
        $work_dir 
    . $species 
    . $opt_type 
    . $contig_opt_flag 
    . ".gp_cds_length";
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
        chomp $line;
        my @columns_f = split " ", $line;
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
    return [
        $temp_gp_cdsgff_contig_eval, 
    $temp_fastagpcontig,
        $temp_tabulargpcontig,       
    $temp_seqlencontig
    ];



}    #processSequences optimization

## GETGENES FUNCTION: EXTRACT FLANKED SEQUENCES FROM GENE MODELS FOR LATER OPTIMIZATION
sub get_genes ($my_fastas_dir, $my_pso_tmp_gp_from_gff, $work_dir)
{
    ## {
    ## BUG last variable is passed empty?
    ## 2016.12.12 my ( $my_fastas_dir, $my_pso_tmp_gp_from_gff, $work_dir ) = @_;

    my $only_non_red = 0;
    my $gp_out_tbl = "$work_dir/gp_exon_utr_intron.tbl";


    open(my $FH_REFGENE, "<", "$my_pso_tmp_gp_from_gff")
        or croak "Failed here";
    ### [<now>] Reading at <file>[<line>]...
    ###  $my_pso_tmp_gp_from_gff
    ## FIXME: $FH_OUT_tblgb to $FH_OUT_tbl_gp
    open(my $FH_OUT_tblgb, ">", "$gp_out_tbl") or croak "Failed here";
    while (<$FH_REFGENE>)
    {

        #my $line = <$FH_REFGENE>;
        #chomp $line;
        my @tmp_cols = split; ## " ";    #, $_;
        ## DEBUG
        #~ print Dump @tmp_cols;
        my (
            $name,     
        $chrom,          
        $strand,
            $tx_start, 
        $tx_end,         
        $cds_start,
            $cds_end,  
        $exon_count_int, 
        $tmp_exons_starts,
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

        my $cds_len    = $cds_end   - $cds_start;
        my $tx_len     = $tx_end    - $tx_start;
        my $cds_offset  = $cds_start - $tx_start;

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
            my @tabular = ();
            while (<$FH_IN>)
            {

                #chomp;
                #print {*STDERR} "$_";
                push @tabular, "$_";
            }
            close $FH_IN;

            #  print {*STDERR} "\nGP: @tabular\n";

            my $sub_seq = "";

            #my $sublen = 0;
            foreach my $line (@tabular)
            {
                chomp $line;
                my @columns_f = split " ", $line;

                #print {*STDERR} "$f[0]\n";
                $sub_seq .= $columns_f[1];

            }
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

            for (my $ii = 0 ; $ii < $exon_count_int ; $ii++)
            ## XXX fixing C style loop XXX
            ##my $ii = 0;
            ##for my $ii (0 .. $exon_count_int)
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

    #$gp_out_tbl = $gp_out_tbl;
    return $gp_out_tbl;

}    #get_genes END

#FIXME: not used. 
#ERROR: Subroutine main::bit_score_graph redefined at (eval 75) line 1 

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


sub BitScoreGraph ($info_output, $info_thresh, $offset)
{
    # 2016.12.13a 
    # {

    # my ($info_output, $info_thresh, $offset) = @_;
    print {*STDERR}
        "bitscoregraph input:  $info_output, $info_thresh, $offset\n";
    my @info = ();    ##($offset - 1, $offset + 1);

    # 2016.12.14a
    # my @fields;

    open(my $FH_INFO, "<", "$info_output") or croak "Failed here";
    while (<$FH_INFO>)
    {    ## 2016.12.14a: these information files do not have anythyng but numbers

        my ($matrix_position, $info_score) = split;
        # @fields = split;
        #~ printf {*STDERR} "%2s %2.2f %s",
        #~ ($fields[0], $fields[1], "=" x int($fields[1] * 30));
        printf {*STDERR} "%2s %2.2f %s",
            ($matrix_position, $info_score, "=" x int($info_score * 30));

        if ($info_score > $info_thresh)
        {
            push(@info, $matrix_position);
        }
        print {*STDERR} "\n";
    }
    close $FH_INFO;

    # DEBUG
    #~ say "INFO SORTING PRE";
    #~ print Dump @info;
    #my @sorted_info = sort { $a <=> $b } @info;
    #~ say "INFO SORTING POST";
    #~ print Dump @sorted_info;

    my $start = (shift @info);
    my $end   = pop @info;
    # 2016.12.15a
    # hack: thresholds may be too strict in some cases???
    # make sure that at least one position bordering the site is included
    if ($start > $offset - 1)
    {
        $start = $offset - 1;
    }
    # 2016.12.15a
    # BUG??
    # it should be one base after the GT/AG/ATG etc 2fix???
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

    print {*STDERR} "\nBitScoreGraph out: $start, $end\n";
    return ($start, $end);
}    #end BitScoreGraph

sub get_pre_matrix ($kmers_tbl, $order)
{  
     
    ## UNUSED
    ## my @orders  = (qw(hmm_0 hmm_2 hmm_3 hmm_4 hmm_5 hmm_6 hmm_7 hmm_8));
    ## my $ordname = $orders[$order];
    my $FH_KMERS;
    open($FH_KMERS, "<", "$kmers_tbl") or croak "Failed here";
    my $first_line = <$FH_KMERS>;
    my ($seq_name, $seq) = split $first_line;
    my $seq_len = length($seq);
    close $FH_KMERS;

    my $frequency_fn = $work_dir . basename($kmers_tbl) . ".freq";
    
    ### [<now>] start frequency.py at <file>[<line>]...
    my $my_command   = ("./bin/frequency.py   1 $kmers_tbl  >  $frequency_fn");
    run($my_command);    
    ### [<now>] finish frequency.py at <file>[<line>]...
}

## GETKMATRIX FUNCTION (Splice sites and Start ATG codon PWMs)
sub get_K_matrix ($true_kmers_tbl, 
          $backgrnd_kmers_tbl, 
          $order, 
          $offset,
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
    
    ## FIXME: not used yet 20200319
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
    ## BUG temp fix
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

    ## BUG?
    my $true_seq_name = $true_kmers_tbl;
    $true_seq_name =~ s/\.tbl$//;
    my $backgrnd_seq_name = $backgrnd_kmers_tbl;
    $backgrnd_seq_name =~ s/\.tbl$//;
    ###  $true_seq_name $backgrnd_seq_name;

    ## Open true sequences
    #    print {*STDERR} "$true_kmers_tbl (true)\n";
    open(my $FH_TRUE_SEQ, "<", "$true_kmers_tbl") or croak "Failed here";
    $_ = <$FH_TRUE_SEQ>;
    my @columns_t = split;
    my $len       = length($columns_t[1]);
    close $FH_TRUE_SEQ;

    ## Open false (background???) sequences
    #    print {*STDERR} "$backgrnd_kmers_tbl (false)\n";
    open(my $FH_BACKGRND_SEQ, "<", "$backgrnd_kmers_tbl")
        or croak "Couldn't open $backgrnd_kmers_tbl: $OS_ERROR \n";
    $_ = <$FH_BACKGRND_SEQ>;
    my @columns_f = split;
    my $len2      = length($columns_f[1]);
    close $FH_BACKGRND_SEQ;

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

    run("./bin/frequency.py   1 $true_kmers_tbl  >  $true_seq_freq_fn");
    run("./bin/frequency.py   1 $backgrnd_kmers_tbl >  $backgrnd_seq_freq_fn");

    my $my_command_A =
        "./bin/information.py    $true_seq_freq_fn $backgrnd_seq_freq_fn ";

    #~ my $my_command_B =
    #~ "| gawk 'NF==2 && \$1<=$my_freq_field_limit_1 && \$1>=$my_freq_field_limit_2'";

    my $my_command_C = " > $my_freq_subtract_fn ";

    #my $my_command = $my_command_A . $my_command_B;
    my $my_command = $my_command_A . $my_command_C;

    say "\n $my_command \n";
    run($my_command);
    ## $temp_infolog = $my_freq_subtract_fn;

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
            "./bin/logratio_zero_order.py   $backgrnd_seq_freq_fn $true_seq_freq_fn > $my_T_generic_logratio_freq_matrix_fn";
        say "\n $my_command \n";
        run($my_command);

    }
    else
    {
        $my_command =
            ##"gawk -f ./bin/Getkmatrix.awk $order $len $true_kmers_tbl | $sort > $my_true_freq_matrix_fn";
            "./bin/get_k_matrix.py   $order $len $true_kmers_tbl | $sort > $my_true_freq_matrix_fn";
        say "\n $my_command \n";
        run($my_command);

        #~ run(
        #~ " gawk -f ./bin/Getkmatrix.awk $order $len $true_kmers_tbl | $sort > $true_seq_name.$ordname-matrix"
        #~ );

        $my_command =
            ##"gawk -f ./bin/Getkmatrix.awk $order $len2 $backgrnd_kmers_tbl | $sort > $my_backgrnd_freq_matrix_fn ";
            "./bin/get_k_matrix.py   $order $len2 $backgrnd_kmers_tbl | $sort > $my_backgrnd_freq_matrix_fn ";
        say "\n $my_command \n";
        run($my_command);

        $my_command =
            ## "gawk -f ./bin/logratio_kmatrix.awk $my_backgrnd_freq_matrix_fn $my_true_freq_matrix_fn > $my_T_generic_logratio_freq_matrix_fn ";
            "./bin/logratio_kmatrix.py   $my_backgrnd_freq_matrix_fn $my_true_freq_matrix_fn > $my_T_generic_logratio_freq_matrix_fn ";

        say "\n $my_command \n";
        run($my_command);

    }

    #need to check output and then go on
    ## draw bit score bar graph function (nested, local)


    my $my_info_thresh = $my_info_thresholds{$matrix_type};

    #say "\n matrix sub: $matrix_type, $my_info_thresh \n";
    my ($start, $end) =
        BitScoreGraph($my_freq_subtract_fn, $my_info_thresh, $offset);
    
    ### $start, $end using $my_freq_subtract_fn \n";

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


    ## DONOR DIMATRIX START
    my $exec_A1 = "./bin/submatrix.py ";
    my $exec_B1 = "./bin/prepare_dimatrix_donor_4parameter.awk ";
    ## my $exec_B1 = "./bin/matrix_4_parameter.py ";

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

        my $my_command_B =
         "$exec_B1 $pre_offset $new_offset $post_offset $my_T_generic_lograt_summatrix_fn > $my_T_generic_matrix_4param_fn";

        # my $my_command_B =
        #    "$exec_B1 $my_T_generic_lograt_summatrix_fn donor $pre_offset $new_offset $post_offset  > $my_T_generic_matrix_4param_fn";

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

        my $my_command_B =
            "./bin/prepare_dimatrix_acceptor_4parameter.awk $pre_offset $new_offset $post_offset $my_T_generic_lograt_summatrix_fn > $my_T_generic_matrix_4param_fn";

        # my $my_command_B =
        #    "$exec_B1 $my_T_generic_lograt_summatrix_fn acceptor $pre_offset $new_offset $post_offset  > $my_T_generic_matrix_4param_fn";

        print {*STDERR} "my_command_B: \n$my_command_B \n";
        run($my_command_B);


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
            "./bin/submatrix.py   $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_lograt_summatrix_fn";
        ##run("gawk -f ./bin/submatrix.awk $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_lograt_summatrix_fn"
        ##);
        say $my_command;
        run($my_command);
        $my_command =
    " ./bin/prepare_trimatrix_start_4parameter.awk $pre_offset $new_offset $post_offset $my_T_generic_lograt_summatrix_fn > $my_T_generic_matrix_4param_fn";
        say $my_command;
        run($my_command);
    
    
        #run(" ./bin/prepare_trimatrix_start_4parameter.awk $pre_offset $new_offset $post_offset $my_T_generic_lograt_summatrix_fn > $my_T_generic_matrix_4param_fn"
        #);
    }

    ## ATG DIMATRIX END

    ## ALL REMAINING CASES START
    else
    {    ### [<now>] Running REMAINING CASES at <file>[<line>]...

        # print {*STDERR} "$path/submatrix_order0.awk $start $end $true_seq_name-log.$ordname-matrix\n";
        $my_command =
            "./bin/submatrix_order0.py  $start $end $my_T_generic_logratio_freq_matrix_fn  > $my_T_generic_matrix_4param_fn";
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

            for (my $i = 0 ; $i < $param->numIsocores ; $i++)
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
                    croak "error in setting exon weights L1915\n";
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
                @evaluation_output = split " ";

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

sub BuildOptimizedParameterFile ( $eval_array )
{
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
    open(my $FH_SOUT, ">", "$work_dir/$species.BuildOptimizedParameterFile.log")
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

    ## FIVE BEST PERFORMANCE RESULTS AFTER OPTIMIZATION
    for (my $i = 0 ; $i < 5 ; $i++)
    {
        print {*STDERR} join("\t", @{$sorted_eval[$i]}), "\n";
    }
    ##

    ## BUILD BEST PERFORMANCE (OPTIMIZED) PARAMETER FILE (USE VALUES OBTAINED ABOVE)

    my $param = Geneid::Param->new();
    $param->readParam("$new_param_fn");

    for (my $i = 0 ; $i < $param->numIsocores ; $i++)
    {
        if (
            !defined @{$param->isocores}[$i]->Exon_weights(
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
            !defined @{$param->isocores}[$i]->Exon_factor(
                    [
                        $best_OlWeight, $best_OlWeight,
                        $best_OlWeight, $best_OlWeight
                    ]
                )
        )
        {
            #     if (!defined @{$param->isocores}[$i]->Exon_factor([0.4,$best_IoWF,$best_IoWF,0.4])) {
            croak "error in setting exon weights\n";
        }
        if (
            !defined @{$param->isocores}[$i]->Site_factor(
                    [
                        1 - $best_OlWeight,
                        1 - $best_OlWeight,
                        1 - $best_OlWeight,
                        1 - $best_OlWeight
                    ]
                )
        )
        {
            #     if (!defined @{$param->isocores}[$i]->Site_factor([0.55,1-$best_IoWF,1-$best_IoWF,0.55])) {
            croak "error in setting exon weights\n";
        }
    }

    #write new parameter file (optimized)
    ## UNUSED
    ## my $optimized_geneid_param_fn = "$results_dir/$species.geneid.optimized.param";
    $param->writeParam("$species.geneid.optimized.param");

    print {*STDERR}
        "\n\nNew optimized parameter file named: $species.geneid.optimized.param \n";
    print {$FH_SOUT}
        "\n\nNew optimized parameter file named: $species.geneid.optimized.param \n";

    close $FH_SOUT;
    return [$best_ExWeightParam, $best_OlWeight, 0, 0, \@evaluation_init];

    ## 2016.12.13c }
    ##close $FH_SOUT;
    ##return 1;
}

sub parameter_evaluate ($gp_fasta, 
                        $gp_gff_fn, 
            $new_param_fn, 
            $OligoWeight_ini,
                        $ExWeightParam_ini)
{
    ### [<now>] Running  parameter_evaluate at <file>[<line>]...
    my $my_command;

    my $geneid_test_predict_gff_fn =
        "$work_dir/geneid_test_predictions." . basename($new_param_fn) . ".gff";
    print {*STDERR}
        "\ngeneid_test_predict_gff_fn: $geneid_test_predict_gff_fn\n";

    my ($FH_GENEID, $fname_geneid) =
        tempfile(DIR => $geneid_dir, SUFFIX => '.eval_par.geneid');
    print "\ntemp geneid file: $fname_geneid \n";
    $my_command = "./bin/geneid -GP $new_param_fn $gp_fasta > $fname_geneid";

    say($my_command);
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
        @evaluation_test = split " ";

        #@evaluation_output = split " ", $_;
    }
    close $FH_IN;

    ### [<now>] Finished sub  parameter_evaluate at <file>[<line>]...


    return \@evaluation_test;

}    # evaluate parameter function

sub calc_stats
{
    ## FIX_ME: 
    ## 20200318
    ## rewrite in the upstream python
    ## output as YAML if needed
    
    ## BUG variable names hard to guess
    my (
        ## 2016.12.15 unused
        ## $species,
        ## $sout,
        $train_introns_filtered_tbl,  $train_cds_filtered_tbl,
        $train_filtered_gff  ,         $train_inframestop_int,
        $eval_inframestop_int,       $train_transcr_used_int,
        $train_noncanon_donors_int,  $train_noncanon_accept_int,
        $train_noncanon_ATGx_int,    $markov_model,
        $total_coding,               $total_noncoding,
        $st_donor,                   $en_donor,
        $st_accept,                  $en_accept,
        $st_ATGx,                    $en_ATGx,

    ) = @_;
    my $my_command = "";

    ## OBTAIN GENE MODEL SET STATISTICS
    ## Open gene model object

    #~ $param->geneModel(Geneid::GeneModel->new());
    #~ $param->geneModel->useDefault;
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
    
    
    ## FIX_ME: move outside 
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
    return ($intron_short_int, $intron_long_int, $intergenic_min,
        $intergenic_max);

}

sub calc_average ($input_numbers)
{
    ##2016.12.13a {
    ## my ($input_numbers) = @_;
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

    return ($mean, $stdev);
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
        my ($seq_name, $seq) = split(/\s+/);
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

    print {*STDERR} "##tbl_2_single_fastas $in_tbl_fn\n";
    while (<$FH_IN_TBL>)
    {
        my $input_line = $_;
        chomp($input_line);

        #my ( $seq_name, $seq ) = split( /\s+/o, $_ );
        my ($seq_name, $seq) = split(/\s+/o, $input_line);

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

    ## 2016.12.12 my ($genetic_code_fn, $cds_fn, $translated_tbl_fn) = @_;
    print {*STDERR} "$genetic_code_fn in loop\n";
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

        my $line = $_;

        my ($aa, $codon) = split(/\s+/, $line);

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
    print {*STDERR} "translating: $cds_fn \n";

    while (<$FH_CDS_IN>)
    {

        my $line = $_;

        my ($name, $seq) = split(/\s+/, $line);

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
    say "\ntranslated_tbl_fn : $translated_tbl_fn \n";
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
            evaluation
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
            prepare_trimatrix_start_4parameter.awk)
    );

    for my $prog_name (@progs_2_check)
    {
        $my_command = "which $prog_name > /dev/null";
        run($my_command);

        #or croak "\n $prog_name not available\n";
    }

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
    $num_seqs, $backgrnd_tbl)
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
            "zcat ./test_data/pmar_1M_60mers.tbl.gz > $backgrnd_tbl";
        run($my_command);
    }

    return 1;

}    # END get_background_kmers function


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
            $my_order = "1";    
            ### [<now>] $my_order Acc/Don site...
        }
    }
    elsif ($my_site_type eq 'ATGx')
    {
        if ($sites_number > $train_sites_markov_cutoff)
        {
            $my_order = "2";    
            ### [<now>] $my_order ATGx site...
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
