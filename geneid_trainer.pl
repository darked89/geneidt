#!/usr/bin/perl

## if run under perlbrew, use i.e.:
#!/$HOME/perl5/perlbrew/perls/perl-5.10.1/bin/perl

## checks & debugs modules
use Modern::Perl;

use strict;
use warnings;
use autodie;
use diagnostics -verbose;
use sigtrap qw(stack-trace old-interface-signals);

use Carp qw(carp cluck croak confess);
use Carp::Always;

#croak "all OK?";

use Carp::Assert qw(assert);

#use Data::Dumper::Perltidy;
use Data::Dumper;

## common modules, used in part at the moment
use Cwd;

#use Env qw(PATH, PERL5LIB);
use Getopt::Long;
use File::Path;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use IPC::System::Simple qw(run system capture EXIT_ANY);
use Readonly;
use feature 'say';

## geneid_trained modules
use Geneid::Param;
use Geneid::Isocore;
use Geneid::geneid;
use Geneid::geneidCEGMA;

## MAIN VARIABLES
my $PROGRAM      = "geneid_trainer";
my $VERSION      = "2016.04.18";
my $PROGRAM_HOME = getcwd;

my $exec_path = "$PROGRAM_HOME/bin/";

local $ENV;
$ENV{'PATH'} = $exec_path . ":" . $ENV{'PATH'};

my $geneticcode = "./etc/genetic.code";

#print STDERR "geneticcode: $geneticcode\n";

## no need to run anything if this fails
check_external_progs();

## Move parts necessary for getting comand line args here
my $species     = "";
my $gff         = "";
my $fasta       = "";
my $sout        = "-";
my $pout        = "-";
my $branchp     = 0;
my $reduced     = 0;
my $interactive = 0;
my $tenfold     = 0;
my $gff2ps      = 0;

my $standard_run = 1;    # get rid of !$reduced

## Get arguments (command line)
GetOptions(
    'species:s'       => \$species,
    'gff:s'           => \$gff,
    'fastas:s'        => \$fasta,
    'sout|statsout:s' => \$sout,

    #'branch'          => \$branchp,
    #'reduced|red'     => \$reduced,

    #'path|binpath:s'  => \$path,
    #'interactive|i' => \$interactive,
    #'tenfold'       => \$tenfold,
    #'gff2ps'        => \$gff2ps
);
my $usage =
"Usage: $0 -species H.sapiens -gff gffname -fastas fastasname -sout <statsfile_out> -branch -reduced -path <executables_path>\n";

print STDERR $usage and exit unless ( $species && $gff && $fasta && $sout );
## EXAMPLE COMMAND LINE: ./geneidTRAINer1_2TA.pl -species S.cerevisiae -gff S_cerevisiae4training.gff -fastas yeast_genome.fa -sout stats.txt -branch -reduced
## Get arguments (command line) END

## end of getting ARGS, do stuff

#~ ## preparing directories
#~ my $work_dir = "$PROGRAM_HOME/00_gtrain_workdir/";
#~ run("mkdir -p $work_dir");
#~ my $tmp_dir = "$work_dir/temp_00/";
#~ run("mkdir -p $tmp_dir;");

#create_work_dirs();

#my $work_dir = "$PROGRAM_HOME/workdir_00_gtrain/";
## NEW changing to relative directory address
my $work_dir = "./workdir_00_gtrain/";

my $tmp_dir = "$work_dir/temp_00/";

#my $TMP         = "$work_dir/tmp/";
#my $stats_dir   = "$work_dir/stats/";
my $stats_dir   = "$work_dir/statistics_${species}/";
my $sites_dir   = "$work_dir/sites/";
my $plots_dir   = "$work_dir/plots/";
my $introns_dir = "$work_dir/introns/";
my $fastas_dir  = "$work_dir/fastas/";
my $backgrd_dir = "$work_dir/backgrd/";
my $cds_dir     = "$work_dir/cds/";
my $geneid_dir  = "$work_dir/geneid/";
my $bigbag_dir  = "$work_dir/bigbag/";

my @data_dirs = (
    $work_dir,  $tmp_dir,    $fastas_dir,  $stats_dir,
    $sites_dir, $plots_dir,  $introns_dir, $backgrd_dir,
    $cds_dir,   $geneid_dir, $bigbag_dir
);

create_data_dirs(@data_dirs);

#my $TMPROOT  = "trainer_$$";
#my $CEGMATMP = "$/tmp_dir$TMPROOT";
my $fh_SOUT;

## PROGRAM SPECIFIC VARIABLES
my $count = 0;

my $temptblcaps = "";
my $label       = 0;

#my $cutoff             = -7;

## Constant values. modify if needed
Readonly::Scalar my $pwm_cutoff     => -7;
Readonly::Scalar my $bases_offset   => 30;    #bases in fron/after? a feature
Readonly::Scalar my $train_fraction => 0.8;   #fraction of seq used for training
Readonly::Scalar my $train_loci_cutoff         => 500;
Readonly::Scalar my $train_sites_cutoff        => 1400;
Readonly::Scalar my $train_sites_cutoff_alt    => 1200;  # changed in some part?
Readonly::Scalar my $train_sites_markov_cutoff => 5500;
## BUG : extracted donor/acceptor/start all are 60bp in size at this point
#Readonly::Scalar my $backgrnd_kmer_size        => 62;
Readonly::Scalar my $backgrnd_kmer_size => 60;
Readonly::Scalar my $backgrnd_kmer_num  => 100000;
## End Constant values
## need to explain or just incorporate Readonly::Scalar
#~ ( $totalcodingbases > 400000 && $totalnoncodingbases > 100000 )
#~ || ( $totalcodingbases > 375000 && $totalnoncodingbases > 150000 )
#~ || (   $totalnoncodingbases > 35000
#~ && $totalcodingbases > ( 25 * $totalnoncodingbases ) )

## another numbers from middle of the code.
#~ if (
#~ !defined @{ $param->isocores }[0]->set_profile(
#~ 'Branch_point_profile', $prof_len_bra, $fxdbraoffset, -50,
#~ $order, 0, 1, 40, 10, 0, 0, $branchmatrix
#~ )

## PROGRAM SPECIFIC VARIABLES (unordered...)

my $pin = "-";

my $answer             = "";
my $validationanswer   = "";
my $jacknifevalidate   = "0";
my $useallseqs         = 0;
my $tempgff4training   = "";
my $tempgff4evaluation = "";
my $templocus_id       = "";
my $templocus_id_new   = "";
my $templocusid_eval   = "";
my $seqsused           = "";
my ( $gpevalgff, $gpevalfa, $gpevaltbl ) = ( "", "", "" );
my ( $gpevalcontiggff, $gpevalcontigfa, $gpevalcontigtbl, $gpevalcontiglen ) =
  ( "", "", "", "" );
my ( $outcdseval, $outintroneval, $outlocus_id_eval, $outgffeval, $inframeeval )
  = ( "", "", "", "", 0 );
my ( $outcds, $outintron, $outlocus_id, $outgff, $inframe ) =
  ( "", "", "", "", 0 );
my $total_seqs      = "";
my $locus_id        = "";
my $new_locus_id    = "";
my $sortn           = "sort -n";
my $seqs4training   = "";
my $seqs4evaluation = "";
my $starttbl        = "";
my $id;
my $tblseq    = "";
my $usebranch = 0;

my $memefile    = "";
my $motifnumber = "";
my $memeanswer  = "";
my ( $branchmatrix, $prof_len_bra, $fxdbraoffset, $startbranch, $endbranch ) =
  ( "", "", "", "", "" );
my @evaluation   = ();
my @jacknifeeval = ();
my ( $bestIeWF, $bestIoWF, $bestAcc, $bestMin ) = ( "", "", "", "" );
my $fullengthbranchtbl = "";

my $geneidgffsorted = "";
my $total_genomic   = "";
my $bckgrnd         = "";

my $reducedtraining = 0;
my (
    $outdonortbl,    $totalnoncandon, $outacceptortbl,
    $totalnoncanacc, $outstarttbl,    $totalnoncansta
) = ( "", "", "", "", "", "" );
my ( $gptraingff, $gptrainfa, $gptraintbl, $gptrainlen, $gpevallen ) =
  ( "", "", "", "", "" );
my ( $gptraincontiggff, $gptraincontigfa, $gptraincontigtbl, $gptraincontiglen )
  = ( "", "", "", "", "" );
my $totalseqs4training = "";
my $templist_train     = "";
my $value              = "";
my $gffseqseval        = "";
my ( $shortintron, $longintron, $minintergenic, $maxintergenic ) =
  ( "", "", "", "" );
my $donorsubprofile    = "";
my $acceptorsubprofile = "";
my $startsubprofile    = "";
my $branchsubprofile   = "";
my $temp_jkf_geneid    = "";

my $contigopt               = 0;
my $tempgeneidgffsorted     = "";
my $tempgeneidgffsortedeval = "";
my $ext_flg                 = 0;
my $bp;
my $optimize = 0;

#############################################################
## INITIAL CHECKS
#############################################################
## TODO
## 1. fasta / gff file accessible?
## 2. limits:
## 2a. >= 500 genes in gff  

## CREATE A VARIABLE SPECIES FOR A GIVEN SPECIES ONLY ONCE####
## XX BUG: it does not store the used variables + subsequent eval $_ is
## XX BUG: against the good PERL practice

# my $varsmemory = $species . ".variables";
# open( my $fh_STORV, ">", "$varsmemory" ) or die;
#############################################################
## Common tasks
#############################################################
## sanity check
#######################################################
## CREATE FASTAS CDS; INTRON, SITES DIRs WITHIN PATH (ONLY FIRST TIME)

#
#######################################################
## store statistics directory variable
my $fh_STORV;
open( $fh_STORV, ">", "$work_dir/geneid_trainer_global.log" )
  or croak "Failed here";

print $fh_STORV Data::Dumper->Dump( [$stats_dir], ['$stats_dir'] );

##CREATE A STATS/PARAMETER FILE
my @timeData = localtime(time);

#STATS DIR CREATED FIRST TIME PIPELINE IS RUN FOR A GIVEN SPECIES
my $statsout = "$stats_dir" . join( '_', @timeData ) . "_$sout";
###OPEN STATISTICS OUTPUT AT THIS TIME...EVERY TIME PIPELINE IS RUN
open( $fh_SOUT, ">", "$statsout" ) or croak "Failed here";
print $fh_SOUT "GENE MODEL STATISTICS FOR $species\n\n";

###########################################################
## CREATING A PARAMETER FILE REGARDLESS OF WHETHER THE TRAINING IS COMPLETE OR SHORT VERSION
########################################
## CREATE BLANK PARAMETER FILE############
my $param = Geneid::Param->new($species);

#set isochores to 1
$param->numIsocores(1);
$param->isocores( [ Geneid::Isocore->new() ] );

## END CREATING A PARAMETER FILE REGARDLESS OF WHETHER THE TRAINING IS COMPLETE OR REDUCED

########################################
## IF THERE IS A MEME BRANCH PROFILE####
########################################

if ($branchp) {
    branch_start();
}

#############################################################
## REDUCED/SHORT TRAINING
## SKIPS ALL BUT BACKGROUND/PWM/MM5
#############################################################
if ($reduced) {
    start_reduced();
}

##########################################
## FULL TRAINING -MANDATORY FOR THE FIRST TIME GENEID IS TRAINED FOR A GIVEN SPECIES
if ( !$reducedtraining ) {
    normal_run();
}

sub normal_run {

## Convert fasta to tabular format
## Fasta process
    print STDERR
      "\nConverting genomics fasta file ($fasta) to tabular format\n";

    my $temptbl = $work_dir . $species . ".genomic.tbl";
    $temptbl = FastaToTbl( $fasta, $temptbl );
    run("sort -o $temptbl $temptbl");

    print STDERR "actg to ACTG conversion of input fasta \n";
    my $tblcaps = "";

    open(
        my $fh_LOCID,
        "-|",
"gawk '{gsub(/_/,\"\",\$1);gsub(/\\./,\"\",\$1);print \$1, toupper(\$2)}' $temptbl "
    );
    while (<$fh_LOCID>) {
        $tblcaps .= $_;
    }
    close $fh_LOCID;

    chomp $tblcaps;
    $temptblcaps = $work_dir . $species . ".genomic.tbl";
    open( my $fh_FOUT_caps, ">", "$temptblcaps" ) or croak "Failed here";
    print $fh_FOUT_caps "$tblcaps";
    close $fh_FOUT_caps;

    #prints a lort of seq on the screen!
    #my @tabular = split(/\n/, $tblcaps);
    #my  @tabular = "";
    #push (@tabular,$tblcaps);
    #foreach my $line (@tabular) {print STDERR "$line END";}

    #    $value = scalar(@tabular) ;
    my $my_command = "gawk '{print \$1}' $temptblcaps | sort | uniq | wc -l";

    my $tot_uniq_fasta_seq = capture($my_command);

    #chomp $value;
    #$value = int($value);

    print STDERR
      "\nThe user has provided $tot_uniq_fasta_seq genomic sequences\n";

    #~ ## store file with number of genomic sequences used in training
    #~ print $fh_STORV Data::Dumper->Dump( [$value], ['$value'] );
    #~ ## store tabular file directory
    #~ print $fh_STORV Data::Dumper->Dump( [$temptbl], ['$temptbl'] );
    #~ ## store CAPPED tabular  file directory
    #~ print $fh_STORV Data::Dumper->Dump( [$temptblcaps], ['$temptblcaps'] );
    #~ ## store fastas dir and plots dir
    #~ print $fh_STORV Data::Dumper->Dump( [$fastas_dir], ['$fastas_dir'] );
    #~ print $fh_STORV Data::Dumper->Dump( [$plots_dir],  ['$plots_dir'] );

## place genomic sequences in "fastas_$species" directory
    print STDERR "move genomic sequences into \"$fastas_dir\" directory\n";
    print STDERR "(also transfer genomic fasta length info)\n\n";
## do not create fastas in diretory if they are already created and their number corresponds to the number of sequences in thr array
## CONVERT GENOMICS FASTA TO MULTI FASTA AND PLACE THEM IN APPROPRIATE DIRECTORY

    print STDERR
"Convert $temptblcaps to multiple genomic fastas and place them in $fastas_dir:\n";

    TblToFastaFile( $fastas_dir, $temptblcaps );
    print STDERR
"\n\nConversion of $temptblcaps to multiple genomic fastas completed..\n\nAdd fasta sequence length information to same directory\n\n";
    write_sizes_from_tbl_fn($temptblcaps);

#################################################

## get locus_id file only first time pipeline is run for a given species #ALL GENE MODELS
    print STDERR "\nEliminate undesirable (_ and .) characters from $gff\n";

    my $filtergff = "";

    $my_command =
"gawk '{OFS=\"\\t\"}{gsub(/\\./,\"\",\$1);gsub(/\\./,\"\",\$9);gsub(/_/,\"\",\$0);print}' $gff";
    $filtergff = capture($my_command);

#~ open $fh_LOCID,
#~ "gawk '{OFS=\"\\t\"}{gsub(/\\./,\"\",\$1);gsub(/\\./,\"\",\$9);gsub(/_/,\"\",\$0);print}' $gff |";
#~ while (<$fh_LOCID>) {
#~ $filtergff .= $_;
#~ }

    ## XXX BUG!!! overwrites input gff!!!
    open( my $fh_FOUT, ">", "$gff" ) or croak "Failed here";
    print $fh_FOUT "$filtergff";
    close $fh_FOUT;

    print STDERR "\nObtain locus_id (list of genomic sequences / genes)\n";

    $my_command = "gawk '{print \$1,\$9}' $gff | sort | uniq ";
    $locus_id   = capture($my_command);

    # say "\n TTT got here TTT\n";
    $templocus_id = $work_dir . $species . "_locus_id";
    open( $fh_FOUT, ">", "$templocus_id" ) or croak "Failed here";
    print $fh_FOUT "$locus_id";
    close $fh_FOUT;
##

## store locus id for a given species the first time we run the pipeline #ALL GENE MODELS
    print $fh_STORV Data::Dumper->Dump( [$templocus_id], ['$templocus_id'] );

## number of gene models TOTAL
    $my_command = " gawk '{print \$2}' $templocus_id | sort | uniq | wc -l";
    $total_seqs = capture($my_command);

#~ $total_seqs =
#~ ` gawk '{print \$2}' $templocus_id | sort | uniq | wc | gawk '{print \$2}' `;
    chomp $total_seqs;
## number of genomic sequences TOTAL
    $my_command    = "gawk '{print \$1}' $templocus_id | sort | uniq | wc -l";
    $total_genomic = capture($my_command);

#~ $total_genomic =
#~ ` gawk '{print \$1}' $templocus_id | sort | uniq | wc | gawk '{print \$1}' `;
    chomp $total_genomic;

    print STDERR
"\nThe gff file ($gff) contains a total of $total_genomic genomic sequences and $total_seqs gene models\n";

## store total number of gene models and genomic sequences containing them
    print $fh_STORV Data::Dumper->Dump(
        [ $total_seqs,   $total_genomic ],
        [ '$total_seqs', '$total_genomic' ]
    );

## get a list of genes TOTAL
    print STDERR "\nObtain list of all genes\n\n";
    my $list_seqs = "";

    $my_command = "gawk '{print \$9}' $gff | sort | uniq ";
    $list_seqs  = capture($my_command);

    my $templist = $work_dir . $species . "_list_all_seqs";
    open( $fh_FOUT, ">", "$templist" ) or croak "Failed here";
    print $fh_FOUT "$list_seqs";
    close $fh_FOUT;

## store list of gene models the first time the pipeline is run for a given species
    print $fh_STORV Data::Dumper->Dump( [$templist], ['$templist'] );

    #XXXVVV

    #RUN ONLY FIRST TIME FOR EACH SPECIES /ONLY FIRST TIME
    if ( $total_seqs >= $train_loci_cutoff || !$tenfold ) {

        $totalseqs4training = int( $train_fraction * $total_seqs );

#print $fh_STORV Data::Dumper->Dump( [$totalseqs4training],             ['$totalseqs4training'] );

        print STDERR
"\nA subset of $totalseqs4training sequences (randomly chosen from the $total_seqs gene models) was used for training\n";
        ## DEBUG shuf => random select
        ## head -$totalseqs4training just the first ones
#my $my_command =           "shuf --head-count=$totalseqs4training $templocus_id | sort | uniq";

        my $my_command =
          "head --lines=$totalseqs4training $templocus_id | sort | uniq";

        $new_locus_id = capture($my_command);

        $templocus_id_new =
          $work_dir . $species . "_locus_id_training_setaside80";
        open( my $fh_FOUT, ">", "$templocus_id_new" ) or croak "Failed here";
        print $fh_FOUT "$new_locus_id";
        close $fh_FOUT;

        print $fh_STORV Data::Dumper->Dump( [$templocus_id_new],
            ['$templocus_id_new'] );

## ASSUMING USER SELECTED TO SET ASIDE SEQUENCES FOR EVALUATION (20%)
        $my_command =
          "gawk '{print \$2}' $templocus_id_new | sort | uniq | wc -l";
        $seqsused = capture($my_command);
        chomp $seqsused;

###################
## gff for training subset
####################
        my $gff4training = "";

        print STDERR
"\nThe new training gff file includes $seqsused gene models (80% of total seqs)\n";
        ## ??? BUG ???
        $my_command =
"gawk '{print \$2\"\$\"}' $templocus_id_new | sort | uniq | egrep -wf - $gff";
        $gff4training = capture($my_command);

        $tempgff4training = $work_dir . $species . ".gff_training_setaside80";
        open( $fh_FOUT, ">", "$tempgff4training" ) or croak "Failed here";
        print $fh_FOUT "$gff4training";
        close $fh_FOUT;

        print $fh_STORV Data::Dumper->Dump(
            [ $tempgff4training,   $seqsused ],
            [ '$tempgff4training', '$seqsused' ]
        );

        print STDERR "\nObtain list of training genes\n\n";

        my $list_seqs_train = "";

        $my_command = "gawk '{print \$9}' $tempgff4training | sort | uniq ";
        $list_seqs_train = capture($my_command);

        $templist_train = $work_dir . $species . "_list_train_seqs_setaside80";
        open( $fh_FOUT, ">", "$templist_train" ) or croak "Failed here";
        print $fh_FOUT "$list_seqs_train";
        close $fh_FOUT;

## Store variable with list of sequences set aside for training
        print $fh_STORV Data::Dumper->Dump( [$templist_train],
            ['$templist_train'] );

#########################
## new locus_id for evaluation test set
#########################
        my $locusideval = "";

        $my_command =
"gawk '{print \$0\"\$\"}' $templist_train | egrep -vwf - $templocus_id";
        $locusideval = capture($my_command);
        chomp $locusideval;

        $templocusid_eval =
          $work_dir . $species . "_locus_id_evaluation_setaside20";
        open( $fh_FOUT, ">", "$templocusid_eval" ) or croak "Failed here";
        print $fh_FOUT "$locusideval";
        close $fh_FOUT;

## Store variable with list of sequences set aside for evaluating
        print $fh_STORV Data::Dumper->Dump( [$templocusid_eval],
            ['$templocusid_eval'] );
#########################
## gff for evaluation test set
#########################

        $my_command =
"gawk '{print \$2\"\$\"}' $templocusid_eval | sort | uniq | egrep -wf - $gff | gawk '{ print \$9}' | sort | uniq | wc -l";
        ## ??? BUG this is a number...
        $gffseqseval = capture($my_command);
        chomp $gffseqseval;

#~ $gffseqseval = ` gawk '{print \$2\"\$\"}' $templocusid_eval | sort | uniq | egrep -wf - $gff | gawk '{ print \$9}' | sort | uniq | wc | gawk '{print \$1}' `;
#~ chomp $gffseqseval;

        print STDERR
"The evaluation gff file includes $gffseqseval gene models (20% of total seqs)\n\n";

        $my_command =
"gawk '{print \$2\"\$\"}' $templocusid_eval | sort | uniq | egrep -wf - $gff ";
        my $gff4evaluation = capture($my_command);

        $tempgff4evaluation =
          $work_dir . $species . ".gff_evaluation_setaside20";
        open( $fh_FOUT, ">", "$tempgff4evaluation" ) or croak "Failed here";
        print $fh_FOUT "$gff4evaluation";
        close $fh_FOUT;

## STORE INFO ON NUMBER OF SEQUENCES TO EVALUATE PLUS GFF FILE  OF SET ASIDE SEQUENCES..
        print $fh_STORV Data::Dumper->Dump(
            [ $tempgff4evaluation,   $gffseqseval ],
            [ '$tempgff4evaluation', '$gffseqseval' ]
        );

    }    # seqs > 500

####LOOP IF WE HAVE FEWER THAN 500 SEQUENCES

    else {    # seqs < $train_loci_cutoff
        ## BUG we do not do jacknife anyway here
        croak "we do not have >= $train_loci_cutoff sequences, quitting now";

        #~ $jacknifevalidate = 1;
        #~ ##ALWAYS USE ALL SEQS IF USER PROVIDES FEWER THAN 500 SEQS
        #~ $useallseqs = 1;

#~ ###STORE INFORMATION AS TO WHETHER WE WOULD LIKE TO RUN JACKNIFE AND THAT SEQS WERE NOT set aside for eval $useallseqs=1 <500 seqs
#~ print $fh_STORV Data::Dumper->Dump(
#~ [ $jacknifevalidate,   $useallseqs ],
#~ [ '$jacknifevalidate', '$useallseqs' ]
#~ );

    }    # seqs < 500

    ######### ONLY EXECUTED FIRST TIME THE PIPELINE IS RUN FOR A GIVEN SPECIES
    #XXXVVV2

    #ONLY FIRST TIME ("NOT SHORT VERSION") FOR A GIVEN SPECIES

    if ( !$useallseqs ) {    ##SET SEQS FOR EVAL AND TRAINING (SUBSETS)

## Convert general gff2 to geneid gff format
## extract and check cds and intron sequences. Remove inframe stops and check all seqs start with ATG and end with STOP
## TRAIN

        print STDERR "\nConvert general gff2 to geneid-gff format\n\n";
        $tempgeneidgffsorted =
          generalGFFtoGFFgeneid( $tempgff4training, $species, ".train" );

        ( $outcds, $outintron, $outlocus_id, $outgff, $inframe ) = @{
            extractCDSINTRON( $tempgeneidgffsorted, $templocus_id_new,
                ".train" )
        };

#  print STDERR " OUTSIDE EXTRACTCDSINTRON outgff: $outgff\noutlocus_id: $outlocus_id\n";
## TRAIN
## EVAL
        $tempgeneidgffsortedeval =
          generalGFFtoGFFgeneid( $tempgff4evaluation, $species, ".eval" );

        (
            $outcdseval, $outintroneval, $outlocus_id_eval,
            $outgffeval, $inframeeval
          )
          = @{
            extractCDSINTRON( $tempgeneidgffsortedeval, $templocusid_eval,
                ".eval" )
          };
###EVAL

    }
    elsif ($useallseqs) {    #USE SAME SEQS TO TRAIN/EVALUATE
## Convert general gff2 to geneid gff format
## extract and check cds and intron sequences. Remove inframe stops and check all seqs start with ATG and end with STOP

        print STDERR
"\nConvert general gff2 to geneid-gff format $gff ###SAME SEQS USED TP TRAIN/EVALUATE\n\n";
        $tempgeneidgffsorted =
          generalGFFtoGFFgeneid( $gff, $species, ".train" );

        ( $outcds, $outintron, $outlocus_id, $outgff, $inframe ) =
          @{ extractCDSINTRON( $tempgeneidgffsorted, $templocus_id, ".train" )
          };

    }    #USE SAME SEQS TO TRAIN/EVALUATE

## extract and check splice sites and start codon. Use only canonical info #IN SEQUENCES USED IN TRAINING
    (
        $outdonortbl,    $totalnoncandon, $outacceptortbl,
        $totalnoncanacc, $outstarttbl,    $totalnoncansta
    ) = @{ extractprocessSITES( $outgff, $outlocus_id ) };

## prepare sequences for optimization of newly developed parameter file (TRAIN)

    print STDERR
"\nConvert gff to gp (golden-path-like)format (artificial contig - concatenated sequences - approx. 800 nt between sequences)\n";

    (
        $gptraincontiggff, $gptraincontigfa,
        $gptraincontigtbl, $gptraincontiglen
    ) = @{ processSequences4Optimization( $outgff, ".train", 1 ) };

    print STDERR
"\nConvert gff to gp (golden-path-like)format (training set for later optimization -400-nt flanked sequences)\n";
    ( $gptraingff, $gptrainfa, $gptraintbl, $gptrainlen ) =
      @{ processSequences4Optimization( $outgff, ".train", 0 ) };
    print STDERR "$gptraingff";

    #~ ## STORE VARIABLE INFO IN DATA DUMPER###
    #~ print $fh_STORV Data::Dumper->Dump(
    #~ [
    #~ $outcds,              $outintron,       $outlocus_id,
    #~ $outgff,              $outcdseval,      $outintroneval,
    #~ $outlocus_id_eval,    $outgffeval,      $inframeeval,
    #~ $tempgeneidgffsorted, $inframe,         $outdonortbl,
    #~ $totalnoncandon,      $outacceptortbl,  $totalnoncanacc,
    #~ $outstarttbl,         $totalnoncansta,  $gptraingff,
    #~ $gptrainfa,           $gptraintbl,      $gptrainlen,
    #~ $gptraincontiggff,    $gptraincontigfa, $gptraincontigtbl,
    #~ $gptraincontiglen
    #~ ],
    #~ [
    #~ '$outcds',           '$outintron',
    #~ '$outlocus_id',      '$outgff',
    #~ '$outcdseval',       '$outintroneval',
    #~ '$outlocus_id_eval', '$outgffeval',
    #~ '$inframeeval',      '$tempgeneidgffsorted',
    #~ '$inframe',          '$outdonortbl',
    #~ '$totalnoncandon',   '$outacceptortbl',
    #~ '$totalnoncanacc',   '$outstarttbl',
    #~ '$totalnoncansta',   '$gptraingff',
    #~ '$gptrainfa',        '$gptraintbl',
    #~ '$gptrainlen',       '$gptraincontiggff',
    #~ '$gptraincontigfa',  '$gptraincontigtbl',
    #~ '$gptraincontiglen'
    #~ ]
    #~ );
########################################

    #NOT USING ALL SEQS FOR TRAINING/EVALUATION ALSO PROCESS EVAL SEQS
    if ( !$useallseqs ) {

## prepare test set for evaluation of newly developed parameter file (EVAL)

        print STDERR
"\nConvert gff to gp (golden-path-like)format (400-nt flanking)(test set for evaluation of new parameter file)\n";
        ( $gpevalgff, $gpevalfa, $gpevaltbl, $gpevallen ) =
          @{ processSequences4Optimization( $outgffeval, ".eval", 0 ) };
        print STDERR "DONE\n";

        print STDERR
"\nConvert gff to gp (golden-path-like)format (test set for evaluation of new parameter file - (artificial contig - concatenated sequences - approx. 800 nt between sequences)\n";
        (
            $gpevalcontiggff, $gpevalcontigfa,
            $gpevalcontigtbl, $gpevalcontiglen
        ) = @{ processSequences4Optimization( $outgffeval, ".eval", 1 ) };
        print STDERR "DONE\n";

## STORE VARIABLE INFO IN DATA DUMPER###
        print $fh_STORV Data::Dumper->Dump(
            [
                $gpevalgff,       $gpevalfa,
                $gpevaltbl,       $gpevallen,
                $gpevalcontiggff, $gpevalcontigfa,
                $gpevalcontigtbl, $gpevalcontiglen,
                $tempgeneidgffsortedeval
            ],
            [
                '$gpevalgff',       '$gpevalfa',
                '$gpevaltbl',       '$gpevallen',
                '$gpevalcontiggff', '$gpevalcontigfa',
                '$gpevalcontigtbl', '$gpevalcontiglen',
                '$tempgeneidgffsortedeval'
            ]
        );
########################################

    }

## DELETE DIRECTORIES NO LONGER NEEDED
    print STDERR
"the CDS-containing directory in $cds_dir is no longer needed..\nRemove directory and its contents\n";

    #rmtree( ["$work_dir/cds/"] );
    print STDERR
"the intron-containing directory in $introns_dir is no longer needed..\nRemove directory and its contents\n";

    #rmtree( ["$work_dir/intron/"] );

#print STDERR "the splice-site-containing directory in $path/sites/ is no longer needed..\nRemove directory and its contents\n";
#	rmtree([ "$path/sites/" ]);
#print STDERR "the fastas-containing directory in $path/fastas_$species/ is no longer needed..\nRemove directory and its contents\n";
#	rmtree([ "$fastas_dir" ]);
###
    close $fh_STORV;

    #XXXVVV2

### EVERYTHING BELOW ALWAYS EXECUTED EVEN ON SHORT VERSION OF THE PIPELINE (REDUCED)

## GET BACKGROUND SEQUENCES

    print
"Obtaining $backgrnd_kmer_num background sequences of $backgrnd_kmer_size nucleotides each for estimating background frequencies of nucleotides\n";

    $bckgrnd = $work_dir . $species . "_background.info";

   #$bckgrnd = getBackground( $kmer, $fasta, $temptblcaps, $numseqs, $bckgrnd );
    getBackground( $backgrnd_kmer_size, $fasta, $temptblcaps,
        $backgrnd_kmer_num, $bckgrnd );
###STORE VARIABLE INFO IN DATA DUMPER###
    #   print STORV Data::Dumper->Dump([$bckgrnd], ['$bckgrnd']);

#############################################################

#########
## get donor site statistics
#########
    my $order = "0";
    my $numbersites;
    $numbersites = num_of_lines_in_file($outdonortbl);

    #my $numbersites =  `wc -l $outdonortbl | gawk '{print \$1}'`;
    #chomp $numbersites;
    #$numbersites = int($numbersites);
    my $donoffset =
      $bases_offset;   #position before intron (last of exon (31) -1 for offset)

    if ( $numbersites > $train_sites_cutoff ) {

        $order = "1";

    }
    elsif ( $numbersites <= $train_sites_cutoff ) {

        $order = "0";
    }

    print STDERR
"\nThere are $numbersites donor sites, enough for a matrix of order $order, prior offset: $donoffset $outdonortbl \n";

    my ( $donormatrix, $prof_len_don, $fxddonoffset, $startdonor, $enddonor ) =
      getKmatrix( $outdonortbl, $bckgrnd, $order, $donoffset, 1, 0, 0, 0, 0, 0,
        0 );
    if (
        !defined @{ $param->isocores }[0]->set_profile(
            'Donor_profile', $prof_len_don, $fxddonoffset, $pwm_cutoff, $order,
            0, 1, 0, 0, 0, 0, $donormatrix
        )
      )
    {
        croak "error in setting profile\n";
    }

    my $donsub = "";

#print STDERR "gawk '{print  substr(\$2,($startdonor-3),($prof_len_don+6))}' $outdonortbl\n";

    $my_command =
"gawk '{print  substr(\$2,($startdonor-3),($prof_len_don+6))}' $outdonortbl ";
    $donsub = capture($my_command);

    $donorsubprofile = $work_dir . $species . ".don.sub.profile";
    open( $fh_FOUT, ">", "$donorsubprofile" ) or croak "Failed here";
    print $fh_FOUT "$donsub";
    close $fh_FOUT;

 #print STDERR "$path/pictogram $donorsubprofile $statsdir/Donor -bits -land\n";
## BUG?
    $my_command =
"./bin/pictogram $donorsubprofile $plots_dir/donor_profile.pictogram -bits -land";
    print "\n$my_command\n";
    run($my_command);

    #unlink $donorsubprofile;

    # print STDERR "donormatrix: \n";
    # foreach my $i (@$donormatrix){
    #     print STDERR join("\t",@$i),"\n";
    #}

#########
## get acceptor site statistics
#########

    $order = "0";

    $numbersites = num_of_lines_in_file($outacceptortbl);

    #$numbersites = `wc -l $outacceptortbl | gawk '{print \$1}'`;
    #chomp $numbersites;
    #$numbersites = int($numbersites);
    print "numbersites in $outacceptortbl: $numbersites\n";
    my $accoffset =
      $bases_offset;   #position after intron (first of exon (31) -1 for offset)

    if ( $numbersites > $train_sites_cutoff ) {

        $order = "1";

    }
    elsif ( $numbersites <= $train_sites_cutoff ) {

        $order = "0";
    }

    print STDERR
"\nThere are $numbersites acceptor sites, enough for a matrix of order $order, offset: $accoffset \n";

    my (
        $acceptormatrix, $prof_len_acc, $fxdaccoffset,
        $startacceptor,  $endacceptor
      )
      = getKmatrix( $outacceptortbl, $bckgrnd, $order, $accoffset, 0, 1, 0, 0,
        0, 0, 0 );
    if (
        !defined @{ $param->isocores }[0]->set_profile(
            'Acceptor_profile', $prof_len_acc, $fxdaccoffset, $pwm_cutoff,
            $order, 0, 1, 0, 0, 0, 0, $acceptormatrix
        )
      )
    {
        croak "error in setting profile\n";
    }

    my $accsub = "";

#print STDERR "gawk '{print  substr(\$2,($startacceptor-3),($prof_len_don+6))}' $outacceptortbl\n";

    $my_command =
"gawk '{print  substr(\$2,($startacceptor-3),($prof_len_acc+6))}' $outacceptortbl ";
    $accsub = capture($my_command);

    $acceptorsubprofile = $work_dir . $species . ".acc.sub.profile";
    open( $fh_FOUT, ">", "$acceptorsubprofile" ) or croak "Failed here";
    print $fh_FOUT "$accsub";
    close $fh_FOUT;

#print STDERR "$path/pictogram $acceptorsubprofile $statsdir/Acceptor -bits -land\n";

## BUG?
    $my_command =
"./bin/pictogram $acceptorsubprofile $plots_dir/acceptor_profile.pictogram -bits -land";
    print "\n$my_command\n";
    run($my_command);

   # run("./bin/pictogram $acceptorsubprofile $plots_dir/Acceptor -bits -land");
   #unlink $acceptorsubprofile;

    #  print STDERR "acceptormatrix: \n";
    #   foreach my $i (@$acceptormatrix){
    #       print STDERR join("\t",@$i),"\n";
    #   }

#########
## get start site statistics
#########

    $order       = "0";
    $numbersites = num_of_lines_in_file($outstarttbl);

    #$numbersites = `wc -l $outstarttbl | gawk '{print \$1}'`;
    #print "#1136 numbersites: $numbersites";
    #chomp $numbersites;
    #$numbersites = int($numbersites);
    my $staoffset =
      $bases_offset;  #before first position of the exon (31)minus 1 for offset)

    if ( $numbersites > $train_sites_markov_cutoff ) {

        $order = "2";

    }
    elsif ( $numbersites <= $train_sites_markov_cutoff ) {

        $order = "0";
    }

    print STDERR
"\nThere are $numbersites start sites, enough for a matrix of order $order, offset: $staoffset \n";

    my ( $startmatrix, $prof_len_sta, $fxdstaoffset, $startstart, $endstart ) =
      getKmatrix( $outstarttbl, $bckgrnd, $order, $staoffset, 0, 0, 1, 0, 0, 0,
        0 );

## write to parameter file
    if (
        !defined @{ $param->isocores }[0]->set_profile(
            'Start_profile', $prof_len_sta, $fxdstaoffset, $pwm_cutoff, $order,
            0, 1, 0, 0, 0, 0, $startmatrix
        )
      )
    {
        croak "error in setting profile\n";
    }
#############################

    my $stasub = "";

    $my_command =
"gawk '{print  substr(\$2,($startstart-3),($prof_len_sta+6))}' $outstarttbl ";
    $stasub = capture($my_command);

    $startsubprofile = $work_dir . $species . ".sta.sub.profile";
    open( $fh_FOUT, ">", "$startsubprofile" ) or croak "Failed here";
    print $fh_FOUT "$stasub";
    close $fh_FOUT;

    $my_command =
"./bin/pictogram $startsubprofile $plots_dir/start_profile.pictogram -bits -land";

    #print "\n$my_command\n";
    run($my_command);

    #unlink $startsubprofile;

    #~ ## OPTIONAL BRANCH STATS (FUNGI NORMALLY, AFTER RUNNING MEME)
    #~ if ($usebranch) {
    #~ refactor_branch_sub();
    #~ }

    #~ ## NO BRANCH AGAIN

## DERIVE INITIAL/TRANSITION MARKOV MODEL

    my ( $markovini, $markovtrans, $totalcoding, $totalnoncoding, $markovmodel )
      = @{ deriveCodingPotential( $outcds, $outintron ) };

    #add markov matrices to the parameter file
    if ( !defined @{ $param->isocores }[0]->Markov_order($markovmodel) ) {
        croak "error in setting Markov_order\n";
    }
    if ( !defined @{ $param->isocores }[0]
        ->Markov_Initial_probability_matrix($markovini) )
    {
        croak "error in setting Markov_Initial_probability_matrix\n";
    }
    if ( !defined @{ $param->isocores }[0]
        ->Markov_Transition_probability_matrix($markovtrans) )
    {
        croak "error in setting Markov_Transition_probability_matrix\n";
    }
######################################

## PRODUCE FILE WITH STATS
    if ($usebranch) {

        ( $shortintron, $longintron, $minintergenic, $maxintergenic ) =
          WriteStatsFile(
            $species,        $sout,           $outintron,
            $outcds,         $outgff,         $inframe,
            $inframeeval,    $seqsused,       $totalnoncandon,
            $totalnoncanacc, $totalnoncansta, $markovmodel,
            $totalcoding,    $totalnoncoding, $startdonor,
            $enddonor,       $startacceptor,  $endacceptor,
            $startstart,     $endstart,       $startbranch,
            $endbranch,      $usebranch,      $useallseqs
          );

    }
    else {

        ( $shortintron, $longintron, $minintergenic, $maxintergenic ) =
          WriteStatsFile(
            $species,        $sout,
            $outintron,      $outcds,
            $outgff,         $inframe,
            $inframeeval,    $seqsused,
            $totalnoncandon, $totalnoncanacc,
            $totalnoncansta, $markovmodel,
            $totalcoding,    $totalnoncoding,
            $startdonor,     $enddonor,
            $startacceptor,  $endacceptor,
            $startstart,     $endstart,
            0,               0,
            0,               $useallseqs
          );

    }

}
print STDERR
"\nshortest intron: $shortintron\nlongest intron: $longintron\nminimum intergenic: $minintergenic\nmaximum intergenic: $maxintergenic\n";

##################################################
## WRITE PRELIMINARY NON-OPTIMIZED PARAMETER FILE
$param->writeParam("$species.geneid.param");
my $newparam = "$species.geneid.param";

##############################################
## Select subset for training/evaluation###
##############################################
#if ( !$reducedtraining )

#~ print STDERR
#~ "\nCHECK: use all sequences ?: $useallseqs\njacknife ?: $jacknifevalidate\n";

##############################################
## CALL SUBS:
##############################################
#if ( !$reducedtraining )

### if reduced training (non-default) do not calculate any of the above ALL OF THE ABOVE MUST BE RUN ONLY FIRST TIME GENEID IS TRAINED FOR A GIVEN SPECIES
###EVERYTHING BELOW WILL BE RUN EVERYTIME THE TRAINING PIPELINE IS RUN WHETHER "REDUCED" OR "FULL"

################################
## OPTIMIZE PARAMETER FILE
################################

print STDERR "\nOptimizing new parameter file\n\n";

## BUG we run non interactive here

my $opttype = "";
if ($interactive) {
    do {
        print STDERR
"\nDo you wish to optimize the internal parameter file on individual 400-nt flanked sequences or on an artificial contig made up of the concatenated flanked sequences (approx. 800 nt between genes)? (Write down \"f(lanked)\" if you choose the first option and \"c(contig)\" if you prefer the second)\n";
        $opttype = readline(STDIN);
    } while ( $opttype !~ /^(flanked|f)|(contig|c)$/i );
}
else {
    $opttype = "contig";
}

if ( $opttype =~ /^(contig|c)$/i ) {

    $contigopt = 1;
    print STDERR
"\nYou chose to optimize the internal parameters of geneid based on the artificial contig and therefore NO 10x cross validation will be performed\n\n";
    my $fh_SOUT;
    open( $fh_SOUT, ">", "$species.middle.log" ) or croak "Failed here";
    print $fh_SOUT
"\nThe user chose to optimize the internal parameters of geneid based on an artificial contig and therefore NO 10x cross validation will be performed\n\n";

    $jacknifevalidate = 0;
}
else {
    $contigopt = 0;
}

## OPTIMIZATION FUNCTION NO BRANCH
my $array_ref = "";

## EXON WEIGHT PARAMETER
my $IeWF = "-4.5";
my $deWF = "0.5";
my $FeWF = "-2.5";
## EXON/OLIGO FACTOR PARAMETER
my $IoWF = "0.25";
my $doWF = "0.05";
my $FoWF = "0.50";
##Minimum Branch Profile Distance
my $iMin = "7";
my $dMin = "2";
my $fMin = "9";
##ACCEPTOR CONTEXT
my $iAccCtx = "40";
my $dAccCtx = "10";
my $fAccCtx = "70";

if ( !$branchp ) {    # no clear separate branch profile
    if ($interactive) {
        my $respo = "";
        do {
            print STDERR
"Use default range values for the optimization of geneid eWF (exon weight) and oWF (exon/oligo factor) internal parameters?\n\n(eWF: $IeWF to $FeWF; step $deWF\noWF: $IoWF to $FoWF; step $doWF)\n\nDo you prefer to change these values? ";
            $respo = readline(STDIN);
        } while ( $respo !~ /^(yes|y)|(n|no)$/i );

        if ( $respo =~ /^(yes|y)/i ) {

            my $sline = "";
            my $eline = "";
            my $dline = "";
            do {
                print STDERR "\nType new initial eWF (IeWF): ";
                $sline = readline(STDIN);
            } while ( $sline !~ /(-*[0-9]*\.*[0-9]+)/ );
            $IeWF = $1;

            do {
                print STDERR "\nType new final eWF (FeWF): ";
                $eline = readline(STDIN);
              } while ( $eline !~ /(-*[0-9]*\.*[0-9]+)/
                || $eline <= $sline );
            $FeWF = $1;

            do {
                print STDERR "\nType step (delta) eWF (deWF)): ";
                $dline = readline(STDIN);
            } while ( $dline !~ /(-*[0-9]*\.*[0-9]+)/ );
            $deWF = $1;

            do {
                print STDERR "\nType new initial oWF (IoWF): ";
                $sline = readline(STDIN);
            } while ( $sline !~ /(-*[0-9]*\.*[0-9]+)/ );
            $IoWF = $1;

            do {
                print STDERR "\nType new final oWF (FoWF): ";
                $eline = readline(STDIN);
              } while ( $eline !~ /(-*[0-9]*\.*[0-9]+)/
                || $eline <= $sline );
            $FoWF = $1;

            do {
                print STDERR "\nType step (delta) oWF (doWF): ";
                $dline = readline(STDIN);
            } while ( $dline !~ /(-*[0-9]*\.*[0-9]+)/ );
            $doWF = $1;

        }
    }
## OPTIMIZATION FUNCTIONS

    if ( !$contigopt ) {
        @evaluation = @{
            OptimizeParameter( $gptrainfa, $gptraingff, $newparam, 0, 0, 0, 0,
                $IeWF, $deWF, $FeWF, $IoWF, $doWF, $FoWF, 0, 0, 0, 0, 0, 0 )
        };

        ( $bestIeWF, $bestIoWF, $bestAcc, $bestMin, $array_ref ) =
          @{ BuildOptimizedParameterFile( \@evaluation, $usebranch, 0, 0, 0 ) };

    }
    elsif ($contigopt) {

        @evaluation = @{
            OptimizeParameter( $gptraincontigfa, $gptraincontiggff,
                $newparam, 0, 0, 0, 0, $IeWF, $deWF, $FeWF, $IoWF, $doWF,
                $FoWF, 0, 0, 0, 0, 0, 0 )
        };

        ( $bestIeWF, $bestIoWF, $bestAcc, $bestMin, $array_ref ) =
          @{ BuildOptimizedParameterFile( \@evaluation, $usebranch, 0, 0, 0 ) };

    }

} ## end of if not branch
elsif ($usebranch) {    #use separate branch profile
    if ($interactive) {
        my $respo = "";
        do {
            print STDERR
"Use automatically selected range values for the optimization of geneid eWF (exon weight)/oWF (exon/oligo factor)/Minimum Branch Distance from Acceptor and Context length in which Branch sites should be scored (AccCtx) internal parameters?\n\n(eWF: $IeWF to $FeWF; step $deWF\noWF: $IoWF to $FoWF; step $doWF\nMinBranchDistance: $iMin to $fMin; step $dMin\nAcceptorBranchCtx: $iAccCtx to $fAccCtx; step $dAccCtx)\n\nDo you prefer to change these values? (we do not recommed you change the Min Branch Distance and AccBranch context)";
            $respo = readline(STDIN);
        } while ( $respo !~ /^(yes|y)|(n|no)$/i );

        if ( $respo =~ /^(yes|y)/i ) {

            my $sline = "";
            my $eline = "";
            my $dline = "";
            do {
                print STDERR "\nType new initial eWF (IeWF): ";
                $sline = readline(STDIN);
            } while ( $sline !~ /(-*[0-9]*\.*[0-9]+)/ );
            $IeWF = $1;

            do {
                print STDERR "\nType new final eWF (FeWF): ";
                $eline = readline(STDIN);
              } while ( $eline !~ /(-*[0-9]*\.*[0-9]+)/
                || $eline <= $sline );
            $FeWF = $1;

            do {
                print STDERR "\nType step (delta) eWF (deWF)): ";
                $dline = readline(STDIN);
            } while ( $dline !~ /(-*[0-9]*\.*[0-9]+)/ );
            $deWF = $1;

            do {
                print STDERR "\nType new initial oWF (IoWF): ";
                $sline = readline(STDIN);
            } while ( $sline !~ /(-*[0-9]*\.*[0-9]+)/ );
            $IoWF = $1;

            do {
                print STDERR "\nType new final oWF (FoWF): ";
                $eline = readline(STDIN);
              } while ( $eline !~ /(-*[0-9]*\.*[0-9]+)/
                || $eline <= $sline );
            $FoWF = $1;

            do {
                print STDERR "\nType step (delta) oWF (doWF): ";
                $dline = readline(STDIN);
            } while ( $dline !~ /(-*[0-9]*\.*[0-9]+)/ );
            $doWF = $1;

            do {
                print STDERR "\nType new initial Min Branch Distance (iMin): ";
                $sline = readline(STDIN);
            } while ( $sline !~ /(-*[0-9]*\.*[0-9]+)/ );
            $iMin = $1;

            do {
                print STDERR "\nType new final Min Branch Distance (fMin): ";
                $eline = readline(STDIN);
              } while ( $eline !~ /(-*[0-9]*\.*[0-9]+)/
                || $eline <= $sline );
            $fMin = $1;

            do {
                print STDERR
                  "\nType step (delta) Min Branch Distance (dMin)): ";
                $dline = readline(STDIN);
            } while ( $dline !~ /(-*[0-9]*\.*[0-9]+)/ );
            $dMin = $1;

            do {
                print STDERR
                  "\nType new initial Acceptor/Branch Context (iAccCtx): ";
                $sline = readline(STDIN);
            } while ( $sline !~ /(-*[0-9]*\.*[0-9]+)/ );
            $iAccCtx = $1;

            do {
                print STDERR
                  "\nType new final Acceptor/Branch Context (fAccCtx): ";
                $eline = readline(STDIN);
              } while ( $eline !~ /(-*[0-9]*\.*[0-9]+)/
                || $eline <= $sline );
            $fAccCtx = $1;

            do {
                print STDERR
                  "\nType step (delta) Acceptor/Branch Context (dAccCtx): ";
                $dline = readline(STDIN);
            } while ( $dline !~ /(-*[0-9]*\.*[0-9]+)/ );
            $dAccCtx = $1;

        }
    }
## OPTIMIZATION FUNCTIONS

    if ( !$contigopt ) {

        @evaluation = @{
            OptimizeParameter(
                $gptrainfa,    $gptraingff,   $newparam,     1,
                $prof_len_bra, $fxdbraoffset, $branchmatrix, $IeWF,
                $deWF,         $FeWF,         $IoWF,         $doWF,
                $FoWF,         $iMin,         $dMin,         $fMin,
                $iAccCtx,      $dAccCtx,      $fAccCtx
            )
        };

        ( $bestIeWF, $bestIoWF, $bestAcc, $bestMin, $array_ref ) = @{
            BuildOptimizedParameterFile(
                \@evaluation,  $branchp, $prof_len_bra,
                $fxdbraoffset, $branchmatrix
            )
        };

    }
    elsif ($contigopt) {

        @evaluation = @{
            OptimizeParameter(
                $gptraincontigfa, $gptraincontiggff, $newparam,
                1,                $prof_len_bra,     $fxdbraoffset,
                $branchmatrix,    $IeWF,             $deWF,
                $FeWF,            $IoWF,             $doWF,
                $FoWF,            $iMin,             $dMin,
                $fMin,            $iAccCtx,          $dAccCtx,
                $fAccCtx
            )
        };

        ( $bestIeWF, $bestIoWF, $bestAcc, $bestMin, $array_ref ) = @{
            BuildOptimizedParameterFile(
                \@evaluation,  $branchp, $prof_len_bra,
                $fxdbraoffset, $branchmatrix
            )
        };

    }

#############################

}

my @evaluationinit = @$array_ref;
my @evaluationtest = ();

############
## EVALUATE PERFORMANCE OF NEW PARAMETER FILE ON TEST SET (IF PRESENT)
############

my $paramopt = "$species.geneid.optimized.param";

if ( !$useallseqs ) {
    my $fh_SOUT;
    open( $fh_SOUT, ">", "$species.useallseqs.log" );

    #print STDERR "CHECK EVALUATE: $gpevalfa, $gpevalgff, $paramopt\n";

    if ( !$contigopt ) {

        @evaluationtest =
          @{ EvaluateParameter( $gpevalfa, $gpevalgff, $paramopt ) };

    }
    elsif ($contigopt) {

        @evaluationtest =
          @{ EvaluateParameter( $gpevalcontigfa, $gpevalcontiggff, $paramopt )
          };
    }

    if ( !$usebranch ) {

        print STDERR
          "\nPerformance of new optimized parameter file on test set:\n\n"
          . join( "\t", @evaluationinit[ 2 .. $#evaluationinit ] ), "\n";

        # print $fh_SOUT
        #     "\nPerformance of new optimized parameter file on test set:\n\n"
        #     . join( "\t", @evaluationinit[ 2 .. $#evaluationinit ] ), "\n";
    }
    elsif ($usebranch) {

        print STDERR
          "\nPerformance of new optimized parameter file on test set:\n\n"
          . join( "\t", @evaluationinit[ 4 .. $#evaluationinit ] ), "\n";

        #my $fh_SOUT;
        #print $fh_SOUT
        #    "\nPerformance of new optimized parameter file on test set:\n\n"
        #    . join( "\t", @evaluationinit[ 4 .. $#evaluationinit ] ), "\n";

    }

    print STDERR join( "\t", @evaluationtest ), "\n\n";

    print $fh_SOUT join( "\t", @evaluationtest ), "\n\n";
    close $fh_SOUT;
}    # if NOT using all seqs for training

#######################################
## END OF MAIN PORTION OF SCRIPT
#######################################

sub extractCDSINTRON {

    my ( $gff, $locus_id, $type ) = @_;

    #ERASE FASTA FILES FOR PARTICULAR SPECIES IF ALREADY EXIST
    #	unlink "$work_dir/cds/${species}.${type}.cds.fa";
    #	unlink "$work_dir/intron/${species}.${type}.intron.fa";

    # #####extract CDS and INTRON SEQUENCES
    #my
    print STDERR "\nEXTRACT CDS and INTRON SEQUENCES from $type set..\n\n";
    open( my $fh_LOCUS, "<", "$locus_id" ) or croak "Failed here";
    print STDERR "$locus_id and $gff\n";
    my $count = 0;
    while (<$fh_LOCUS>) {
        my ( $genomic_id, $gene_id ) = split;
        run(" egrep -w '$gene_id\$' $gff > $tmp_dir/$gene_id.gff");
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
    my $tempcdsfa = $cds_dir . ${species} . "$type" . ".cds.fa";
    print STDERR "$tempcdsfa\n\n";
    my $tempcds = $cds_dir . ${species} . "$type" . ".cds.tbl";
    $tempcds = FastaToTbl( $tempcdsfa, $tempcds );
    print STDERR "cds tabular file created for $type sequences \n";

    # ##INTRON
    my $tempintronfa = $introns_dir . ${species} . "$type" . ".intron.fa";
    my $tempintron   = $introns_dir . ${species} . "$type" . ".intron.tbl";
    $tempintron = FastaToTbl( $tempintronfa, $tempintron );

## INTRONS LARGER THAN 0 ONLY

    my $introntblpositive = "";
    my $my_command = "gawk '{if(length(\$2)>0){print \$1,\$2}}' $tempintron ";
    $introntblpositive = capture($my_command);

    #~ open( my $fh_LOCID,
    #~ "gawk '{if(length(\$2)>0){print \$1,\$2}}' $tempintron |" );
    #~ while (<$fh_LOCID>) {

    #~ $introntblpositive .= $_;
    #~ }
    #~ close $fh_LOCID;

    my $tempallintron_positive =
      $work_dir . $species . "$type" . ".intron_positivelength.tbl";

    open( my $fh_FOUT, ">", "$tempallintron_positive" ) or croak "Failed here";
    print $fh_FOUT "$introntblpositive";
    close $fh_FOUT;

    print STDERR
      "intron tabular file created with introns with more than 0 nucleotides\n";
## GET LIST OF SEQUENCES WITH LENGTH >0 and EXCLUDE FROM CDS/locus_id/gff FILES SEQUENCES WITH INTRONS WITH 0 LENGTH
    my $intronzero = "";
    $my_command =
"gawk '{if(length(\$2)==0){print \$1}}' $tempintron | sed 's/\\(.*\\)\\..*/\\1\\_/' | sort | uniq ";
    $intronzero = capture($my_command);

    my $tempall_intron_zero_list =
      $work_dir . $species . "$type" . ".intron_zerolength.list";

    open( $fh_FOUT, ">", "$tempall_intron_zero_list" );
    print $fh_FOUT "$intronzero";
    close $fh_FOUT;

    my $intronzero2 = "";
    $my_command =
"gawk '{if(length(\$2)==0){print \$1}}' $tempintron | sed 's/\\(.*\\)\\..*/\\1/' | sort | uniq ";
    $intronzero2 = capture($my_command);

    my $tempall_intron_zero_list2 =
      $work_dir . $species . "$type" . ".intron_zerolength.list2";

    open( $fh_FOUT, ">", "$tempall_intron_zero_list2" );
    print $fh_FOUT "$intronzero2";
    close $fh_FOUT;

## FILTER SEQUENCES WITH 0 SIZE INTRONS FROM CDS!

    $my_command = "egrep -vf $tempall_intron_zero_list $tempcds ";
    my $cdstblnozero = capture($my_command);

    my $tempallcds_nozero = $work_dir . $species . "$type" . ".cds_nozero.tbl";

    open( $fh_FOUT, ">", "$tempallcds_nozero" );
    print $fh_FOUT "$cdstblnozero";
    close $fh_FOUT;
## ENSURE LOCUSID DOES NOT CONTAIN SEQUENCES WITH 0 SIZE INTRONS
    $my_command = "egrep -vwf $tempall_intron_zero_list2 $locus_id ";
    my $locusidnozero = capture($my_command);

    my $templocus_id_nozero =
      $work_dir . $species . "$type" . "_locus_id_nozero";
    open( $fh_FOUT, ">", "$templocus_id_nozero" );
    print $fh_FOUT "$locusidnozero";
    close $fh_FOUT;
## ENSURE GFF DOES NOT CONTAIN SEQUENCES WITH 0 SIZE INTRONS

    my $gffnozero = "";
    $my_command = "egrep -vwf $tempall_intron_zero_list2 $gff ";
    $gffnozero  = capture($my_command);

    my $tempgffnozero = $work_dir . $species . "$type" . ".nozero.gff";

    open( $fh_FOUT, ">", "$tempgffnozero" ) or croak "Failed here";
    print $fh_FOUT "$gffnozero";
    close $fh_FOUT;

    #	rmtree([ "$path/cds/" ]);
    #	rmtree([ "$path/intron/" ]);

## Convert sequences to protein format and check for in-frame stops
    print STDERR
"\nConvert sequences to protein format and check for in-frame stops and for proteins not starting with an M or not ending with a STOP\n\n";

## SHOWS WHERE GENETIC CODE FILE IS LOCATED AND ITS NAME

    my $tempall_protein = $work_dir . $species . "$type" . ".protein";

    # $tempall_protein = Translate($geneticcode,$tempcds,$tempall_protein);
    $tempall_protein =
      Translate( $geneticcode, $tempallcds_nozero, $tempall_protein );

    $my_command =
"gawk '{print \$2,\$1}' $tempall_protein | egrep '[A-Z]\\*[A-Z]\|^[^M]\|[^\\*] ' | gawk '{print \$2}' | wc -l";
    ## BUG reports just the number
    my $inframestops = capture($my_command);
    chomp $inframestops;

#~ #my $inframestops =`gawk '{print \$2,\$1}' $tempall_protein | egrep '[A-Z]\\*[A-Z]\|^[^M]\|[^\\*] ' | gawk '{print \$2}' | wc | gawk '{print \$1}'`;
#~ chomp $inframestops;

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
        print $fh_FOUT "$inframe";
        close $fh_FOUT;
## REMOVE SEQUENCES WITH IN-FRAME STOPS FROM ORIGINAL CDS / INTRON / LOCUS_ID /GFF FILES AND PRINT NEW FILES
        print STDERR
"\nremove sequences with in-frame stop signals from cds/intron files\n\n";

        $my_command =
"sed 's/\\(.*\\)/\\1_/g' $inframe_protein | egrep -vf - $tempallcds_nozero";
        my $cdstbl2 = capture($my_command);

        my $tempall_cds2 = $work_dir . $species . "$type" . ".cds_filter1.tbl";

        open( $fh_FOUT, ">", "$tempall_cds2" );
        print $fh_FOUT "$cdstbl2";
        close $fh_FOUT;

        my $introntbl2 = "";

        $my_command =
"sed 's/\\(.*\\)/\\1\.i/g' $inframe_protein | egrep -vf - $tempallintron_positive ";
        $introntbl2 = capture($my_command);

        my $tempall_intron2 =
          $work_dir . $species . "$type" . ".intron_filter1.tbl";

        open( $fh_FOUT, ">", "$tempall_intron2" );
        print $fh_FOUT "$introntbl2";
        close $fh_FOUT;

        $my_command =
"sed 's/\\(.*\\)/\\1\$/g' $inframe_protein | egrep -vf - $templocus_id_nozero ";
        my $new_locus_id_filter1 = capture($my_command);

        my $templocus_id_new2 =
          $work_dir . $species . "$type" . "_locus_id_filter_noinframe";

        open( $fh_FOUT, ">", "$templocus_id_new2" );
        print $fh_FOUT "$new_locus_id_filter1";
        close $fh_FOUT;

        #my $gffnew = "";
        $my_command =
"sed 's/\\(.*\\)_.*/\\1\$/g' $inframe_protein | egrep -vf - $tempgffnozero ";
        my $gffnew = capture($my_command);

        my $tempnewgff = $work_dir . $species . "$type" . ".noinframe.gff";

        open( $fh_FOUT, ">", "$tempnewgff" ) or croak "Failed here";
        print $fh_FOUT "$gffnew";
        close $fh_FOUT;

######
        return [
            $tempall_cds2, $tempall_intron2, $templocus_id_new2,
            $tempnewgff,   $inframestops
        ];

    }
    else {    ## ??? END IF THERE ARE INFRAME STOPS
        return [
            $tempallcds_nozero,   $tempallintron_positive,
            $templocus_id_nozero, $tempgffnozero,
            0
        ];
    }    #######END ELSE IF NO SEQS  ARE INFRAME

}    #########sub extractCDSINTRON

## FUNCTION TO EXTRACT AND PROCESS SPLICE SITES AND START CODON
sub extractprocessSITES {

    my ( $gff, $locus_id ) = @_;
    my $fh_LOCID;

## SPLICE SITES
    print STDERR "\nEXTRACT START AND SPLICE SITES from transcripts\n\n";

    #print STDERR "$locus_id and $gff\n";
    my @newsites = ();
    my $count    = 0;

    open( my $fh_LOC_sites, "<", "$locus_id" ) or croak "Failed here";
    while (<$fh_LOC_sites>) {
        my ( $genomic_id, $gene_id ) = split;

        #  print STDERR "$genomic_id,$gene_id\n";
        run("egrep -w '$gene_id\$' $gff > $tmp_dir/$gene_id.gff");

        #  print STDERR "$gene_id $gff $tmp_dir/$gene_id.gff \n\n";
        ## POTENTIAL BUG SPLIT
        my $my_command =
"./bin/ssgff -dabeE $fastas_dir/$genomic_id $tmp_dir/$gene_id.gff > $tmp_dir/${gene_id}.all_sites";
        run($my_command);

#   ` ./bin/ssgff -dabeE $fastas_dir/$genomic_id $tmp_dir/$gene_id.gff > $tmp_dir/${gene_id}.all_sites`;
        foreach my $site (qw(Acceptor Donor Stop Start)) {

#	print STDERR "egrep -A 1 $site $tmp_dir/${gene_id}.all_sites $sitesdir/${site}_sites.fa\n";
## POTENTIAL BUG, split command below

            run(
" egrep -A 1 $site $tmp_dir/${gene_id}.all_sites | sed -e '/--/d' -e '/^\$/d' >> $sites_dir/${site}_sites.fa"
            );
        }
        $count++;
        print STDERR "$count..";
    }    #while $fh_LOC_sites
    close $fh_LOC_sites;

    my $accsites = "$sites_dir/Acceptor_sites.fa";

    #print STDERR "$accsites\n..";
    my $donsites   = "$sites_dir/Donor_sites.fa";
    my $startsites = "$sites_dir/Start_sites.fa";
    my $stopsites  = "$sites_dir/Stop_sites.fa";

    my $prestarttbl = "$work_dir/Start_sites.tbl";
    $prestarttbl = FastaToTbl( $startsites, $prestarttbl );
    my $acceptortbl = "$work_dir/Acceptor_sites.tbl";
    $acceptortbl = FastaToTbl( $accsites, $acceptortbl );

    #print STDERR "$acceptortbl\n";
    my $donortbl = "$work_dir/Donor_sites.tbl";
    $donortbl = FastaToTbl( $donsites, $donortbl );

##ADD N TO START SITES############
    ## POTENTIAL BUG
    my $my_command =
"gawk '{printf \$1\" \";for (i=1;i<=60-length(\$2);i++) printf \"n\"; print \$2}' $prestarttbl > $sites_dir/Start_sites_complete.tbl";
    run($my_command);

#`gawk '{printf \$1" ";for (i=1;i<=60-length(\$2);i++) printf "n"; print \$2}' $prestarttbl > $sites_dir/Start_sites_complete.tbl`;
    my $starttbl = "$sites_dir" . "Start_sites_complete.tbl";
#################################

    print STDERR "\n\nEliminate non-canonical donors/acceptors/starts:\n";

    #  	##EXTRACT NON CANONICAL DONORS
    #my $noncanonical     = "";
    my $noncanonicalname = "";
    my $totnoncanonical  = "";
    my $totcanonical     = "";
    my $newdonortbl      = "";

    $my_command = "gawk '{print \$2}' $donortbl  | egrep -v '^[NATCGn]{31}GT' ";
    my $noncanonical = capture($my_command);

    my $tempdonornoncanonical = $work_dir . $species . "_non_canonical_donor";
    open( my $fh_FOUT, ">", "$tempdonornoncanonical" ) or croak "Failed here";
    print $fh_FOUT "$noncanonical";
    close $fh_FOUT;

    $totnoncanonical = num_of_lines_in_file($tempdonornoncanonical);

    print STDERR
"\nThere are $totnoncanonical non-canonical donors within the training set:\n";

###########################
    if ($totnoncanonical) {    #if there are non canonical donors

        my @noncanonicalname = ();
        open $fh_LOCID,
"egrep -wf $tempdonornoncanonical $donortbl | gawk '{print \$1}' - | sort | uniq |";
        while (<$fh_LOCID>) {

            push( @noncanonicalname, "$_" );

        }
        close $fh_LOCID;

        foreach my $line (@noncanonicalname) {

            #   my (@noncan)= split (/\.\d+:/, $line);
            my (@noncan) = split( /:/, $line );
            my $first = $noncan[0] . ":";
            $noncanonicalname .= "$first\n";

        }

        #  unlink $tempdonornoncanonical;

        my $tempnoncanonicalname =
          $work_dir . $species . "_non_canonical_donor_seq_name";
        open( my $fh_FOUT, ">", "$tempnoncanonicalname" )
          or croak "Failed here";
        print $fh_FOUT "$noncanonicalname";
        close $fh_FOUT;

        open $fh_LOCID, "egrep -vf $tempnoncanonicalname $donortbl |";
        while (<$fh_LOCID>) {
            $newdonortbl .= $_;
        }
        close $fh_LOCID;

        my $tempcanonicaldonor = $work_dir . $species . ".canonical.donor.tbl";
        open( $fh_FOUT, ">", "$tempcanonicaldonor" ) or croak "Failed here";
        print $fh_FOUT "$newdonortbl";
        close $fh_FOUT;

        # unlink $tempnoncanonicalname;

        $totcanonical = num_of_lines_in_file($tempcanonicaldonor);

        print STDERR
"\nThere are $totcanonical canonical donors within the training set:\n";

        push( @newsites, "$tempcanonicaldonor" );
        push( @newsites, "$totnoncanonical" );
    }
    else {    #if there are no non-canonical
        my $totcanonical = num_of_lines_in_file($donortbl);

        print STDERR
          "There are $totcanonical canonical donors within the training set:\n";
        push( @newsites, "$donortbl" );
        push( @newsites, "" );

    }    #if there are no non-canonical

    #  	###########################

    #  	####
    #  	##EXTRACT NON CANONICAL ACCEPTORS
    $noncanonical     = "";
    $noncanonicalname = "";
    $totcanonical     = "";
    my $newacceptortbl = "";
    my $foobar_tmp     = "";

    $my_command =
      "gawk '{print \$2}' $acceptortbl | egrep -v '^[NATCGn]{28}AG'";

    # BUG this blows if there are no such sites...
    #$foobar_tmp = capture($my_command);

    open $fh_LOCID,
      "gawk '{print \$2}' $acceptortbl | egrep -v '^[NATCGn]{28}AG' |";
    while (<$fh_LOCID>) {
        $noncanonical .= $_;
    }
    if ( length($noncanonical) > 0 ) {
        close $fh_LOCID;
    }

    my $tempacceptornoncanonical =
      $work_dir . $species . "_non_canonical_acceptor";
    open( $fh_FOUT, ">", "$tempacceptornoncanonical" ) or croak "Failed here";
    print $fh_FOUT "$noncanonical";
    close $fh_FOUT;

    $totnoncanonical = num_of_lines_in_file($tempacceptornoncanonical);

    print STDERR
"\nThere are $totnoncanonical non-canonical acceptors within the training set:\n";
###########################
    if ($totnoncanonical) {    #if there are non-canonical acceptors

        my @noncanonicalname = ();
        open $fh_LOCID,
"egrep -f $tempacceptornoncanonical $acceptortbl | gawk '{print \$1}' - | sort | uniq |";
        while (<$fh_LOCID>) {
            push( @noncanonicalname, "$_" );
        }

        close $fh_LOCID;

        foreach my $line (@noncanonicalname) {
            my (@noncan) = split( /:/, $line );
            my $first = $noncan[0] . ":";
            $noncanonicalname .= "$first\n";

        }

        unlink $tempacceptornoncanonical;

        my $tempnoncanonicalname =
          $work_dir . $species . "_non_canonical_acceptor_seq_name";

        open( my $fh_FOUT, ">", "$tempnoncanonicalname" );
        print $fh_FOUT "$noncanonicalname";
        close $fh_FOUT;

        open $fh_LOCID, "egrep -vf $tempnoncanonicalname $acceptortbl |";
        while (<$fh_LOCID>) {
            $newacceptortbl .= $_;
        }
        close $fh_LOCID;

        #unlink $tempnoncanonicalname;

        my $tempcanonicalacceptor =
          $work_dir . $species . ".canonical.acceptor.tbl";
        open( $fh_FOUT, ">", "$tempcanonicalacceptor" ) or croak "Failed here";
        print $fh_FOUT "$newacceptortbl";
        close $fh_FOUT;

        #unlink $tempnoncanonicalname;

        my $totcanonical = num_of_lines_in_file($tempcanonicalacceptor);

        print STDERR
"\nThere are $totcanonical canonical acceptors within the training set:\n";

        push( @newsites, "$tempcanonicalacceptor" );
        push( @newsites, "$totnoncanonical" );

    }
    else {    #if there are only canonical use initial file list
        my $totcanonical = num_of_lines_in_file($acceptortbl);

        print STDERR
"There are $totcanonical canonical acceptors within the training set:\n";
        push( @newsites, "$acceptortbl" );
        push( @newsites, "0" );
    }    #if there are only canonical use initial file list

    #  	###########################

    #  	###
    #  	##EXTRACT NON CANONICAL STARTS

    $noncanonical     = "";
    $noncanonicalname = "";
    $totcanonical     = "";
    my $newstarttbl = "";

    open $fh_LOCID,
      "gawk '{print \$2}' $starttbl | egrep -v '^[NATCGn]{30}ATG' |";
    while (<$fh_LOCID>) {
        $noncanonical .= $_;
    }
    if ( length($noncanonical) > 0 ) {
        close $fh_LOCID;
    }

    #close $fh_LOCID;

    my $tempstartnoncanonical = $work_dir . $species . "_non_canonical_start";
    open( $fh_FOUT, ">", "$tempstartnoncanonical" ) or croak "Failed here";
    print $fh_FOUT "$noncanonical";
    close $fh_FOUT;

    $totnoncanonical = num_of_lines_in_file($tempstartnoncanonical);

    print STDERR
"\nThere are $totnoncanonical non-canonical starts within the training set:\n";
###########################

    if ($totnoncanonical) {    #if there are non-canonical starts

        my @noncanonicalname = ();
        open $fh_LOCID,
"egrep -wf $tempstartnoncanonical $starttbl | gawk '{print \$1}' - | sort | uniq |";
        while (<$fh_LOCID>) {
            push( @noncanonicalname, "$_" );
        }
        close $fh_LOCID;

        foreach my $line (@noncanonicalname) {
            my (@noncan) = split( /:/, $line );
            my $first = $noncan[0] . ":";
            $noncanonicalname .= "$first\n";

        }

        unlink $tempstartnoncanonical;

        my $tempnoncanonicalname =
          $work_dir . $species . "_non_canonical_start_seq_name";
        open( $fh_FOUT, ">", "$tempnoncanonicalname" );
        print $fh_FOUT "$noncanonicalname";
        close $fh_FOUT;

        open( $fh_LOCID, "egrep -vf $tempnoncanonicalname $starttbl |" );
        while (<$fh_LOCID>) {
            $newstarttbl .= $_;
        }
        close $fh_LOCID;

        # unlink $tempnoncanonicalname;

        my $tempcanonicalstart = $work_dir . $species . ".canonical.start.tbl";
        open( $fh_FOUT, ">", "$tempcanonicalstart" ) or croak "Failed here";
        print $fh_FOUT "$newstarttbl";
        close $fh_FOUT;

        #unlink $tempnoncanonicalname;
        my $totcanonical = num_of_lines_in_file($tempcanonicalstart);

        print STDERR
"\nThere are $totcanonical canonical starts within the training set:\n";

        push( @newsites, "$tempcanonicalstart" );
        push( @newsites, "$totnoncanonical" );

    }
    else {
        my $totcanonical = num_of_lines_in_file($starttbl);

        print STDERR
"\nThere are $totcanonical canonical starts within the training set:\n";
        push( @newsites, "$starttbl" );
        push( @newsites, "0" );

    }    #if there are only canonical starts
###########################

    return \@newsites;

}    #subectractprocesssites

## FUNCTION TO OBTAIN MARKOV MODELS CORRESPONDING TO THE CODING POTENTIAL
sub deriveCodingPotential {
    my ( $cds, $intron ) = @_;

    my $markov  = "";
    my $markovm = "";

    my $my_command       = "gawk '{ l=length(\$2); L+=l;} END{ print L;}' $cds";
    my $totalcodingbases = capture($my_command);

    #~ ` gawk '{ l=length(\$2); L+=l;} END{ print L;}' $cds `;
    chomp $totalcodingbases;

    $my_command = "gawk '{ l=length(\$2); L+=l;} END{ print L;}' $intron ";
    my $totalnoncodingbases = capture($my_command);

    #~ ` gawk '{ l=length(\$2); L+=l;} END{ print L;}' $intron `;
    chomp $totalnoncodingbases;

    print STDERR
"There are $totalcodingbases coding bases and $totalnoncodingbases non-coding bases on this training set:\n";

    if (
           ( $totalcodingbases > 400000 && $totalnoncodingbases > 100000 )
        || ( $totalcodingbases > 375000 && $totalnoncodingbases > 150000 )
        || (   $totalnoncodingbases > 35000
            && $totalcodingbases > ( 25 * $totalnoncodingbases ) )
      )
    {
        $markov  = "5";
        $markovm = "4";
        print STDERR "Deriving a markov model of order $markov OPTION_1\n";

    }
    else {
        $markov  = "5";
        $markovm = "4";
        print STDERR "Deriving a markov model of order $markov  OPTION_2\n";
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

    print STDERR "Intron model\n markov: ($markovm)";

    my $intron_initial =
      geneidCEGMA::SequenceModel->new( 'intron', 'FREQ', $markovm,
        \@introndois, 10, 0 );

#my $intron_initial = new geneidCEGMA::SequenceModel('intron', 'FREQ', $markovm, \@introndois, 10, 0 );

    #print STDERR "($markovm) - intron initial: $intron_initial\n";
    # $intron_initial->write("$DIR/intron.initial.5.freq");

    my $intron_transition =
      geneidCEGMA::SequenceModel->new( 'intron', 'MM', $markov, \@introndois,
        10, 0 );

#my $intron_transition = new geneidCEGMA::SequenceModel( 'intron', 'MM', $markov, \@introndois, 10, 0 );

    # $intron_transition->write("$DIR/intron.transition.5.freq");

    print STDERR "Coding model\n";

    my $coding_initial =
      geneidCEGMA::SequenceModel->new( 'coding', 'FREQ', $markov - 1,
        \@coding, 0.25, 2 );

#my $coding_initial = new geneidCEGMA::SequenceModel( 'coding', 'FREQ', $markov - 1, \@coding, 0.25, 2 );

    #$coding_initial->write("$path/coding.initial.5.freq");

    my $coding_transition =
      geneidCEGMA::SequenceModel->new( 'coding', 'MM', $markov, \@coding, 0.25,
        2 );

#my $coding_transition = new geneidCEGMA::SequenceModel( 'coding', 'MM', $markov, \@coding, 0.25, 2 );

    #$coding_transition->write("$path/coding.transition.5.freq");

    my $initial_logs =
      geneidCEGMA::log_ratio( $coding_initial, $intron_initial );

    my $transition_logs =
      geneidCEGMA::log_ratio( $coding_transition, $intron_transition );

    geneidCEGMA::write_log( $initial_logs, "$work_dir/coding.initial.5.logs" );

    geneidCEGMA::write_log( $transition_logs,
        "$work_dir/coding.transition.5.logs" );

    #  print STDERR "\nINITIAL LOGS:".${$initial_logs}."\n";
    #open (PROF1,"<$tempcdsintroninitgeneid");
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

    #open (PROF2,"<$tempcdsintrontrangeneid");
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
        \@profileinit,        \@profiletran, $totalcodingbases,
        $totalnoncodingbases, $markov
    ];

}    #derive coding potential

## PROCESS SEQUENCES FUNCTION ( FLANKED GENE MODELS OBTAINED FOR OPTIMIZATION)
sub processSequences4Optimization {

    my ( $gff, $type, $contigopt ) = @_;

    my $outtblname = "";
    my $tblgp      = "";
    my $gff2gp     = "";
    my $fastagp    = "";
    my $gffgp      = "";
    my $my_command = "";

    #my $work_dir;

    open( my $fh_LOCID, "./bin/gff2gp.awk $gff | sort -k 1 |" );
    while (<$fh_LOCID>) {

        $gff2gp .= $_;
    }
    close $fh_LOCID;

    my $tempgff2gp = $work_dir . $species . $type . ".gp";

    open( my $fh_FOUT, ">", "$tempgff2gp" );
    print $fh_FOUT "$gff2gp";
    close $fh_FOUT;
    print STDERR
      "BEFORE GETGENES: $fastas_dir,$tempgff2gp,$work_dir/,$outtblname\n";
    my $pretblgp = GetGenes( $fastas_dir, $tempgff2gp, $work_dir, $outtblname );
    print STDERR "PRETBL AFTER GETGENES: $pretblgp \n";

    print STDERR
"\nGet sequences of 400-nt flanked sequences in tabular and gff formats\n";

#~ my $seq4Optimization_temp_1_fn =  "$work_dir/processSequences4Optimization_temp1.txt";
#~ $my_command =  "gawk 'BEGIN{b=\"x\"}{if (\$1!=b){printf \"\\n\%s \",\$1}printf \"\%s\",\$6;b=\$1}END{printf \"\\n\"}' $pretblgp > $seq4Optimization_temp_1_fn";
#~ run($my_command);
#~ $my_command =  "sort seq4Optimization_temp_1_fn | sed '/^\$/d' - | gawk '{print \$1, toupper(\$2)}' - |";
#open( $fh_LOCID, $my_command );
    open( $fh_LOCID,
"gawk 'BEGIN{b=\"x\"}{if (\$1!=b){printf \"\\n\%s \",\$1}printf \"\%s\",\$6;b=\$1}END{printf \"\\n\"}' $pretblgp | sort | sed '/^\$/d' - | gawk '{print \$1, toupper(\$2)}' - |"
    );

    while (<$fh_LOCID>) {
        $tblgp .= $_;
    }
    close $fh_LOCID;

    my $tempgp_tbl = $work_dir . $species . $type . ".gp.tbl";
    open( $fh_FOUT, ">", "$tempgp_tbl" ) or croak "Failed here";
    print $fh_FOUT "$tblgp";
    close $fh_FOUT;

    open( $fh_LOCID,
"gawk 'BEGIN{OFS=\"\\t\";pos=1;b=\"x\"}{if (\$1!=b){pos=1}; print \$1,\"annotations\",\$3,pos,pos+\$5-1,\"\.\",\"+\",\"\.\",\$1\$2; pos+=\$5;b=\$1 }' $pretblgp | egrep -v '(Intron|Utr)' - |"
    );
    while (<$fh_LOCID>) {
        $gffgp .= $_;
    }
    close $fh_LOCID;

    my $tempgp_gff = $work_dir . $species . $type . ".gp.gff";
    open( $fh_FOUT, ">", "$tempgp_gff" ) or croak "Failed here";
    print $fh_FOUT "$gffgp";
    close $fh_FOUT;

    print STDERR "DONE\n";

    print STDERR
      "\nGet sequences of 400-nt flanked sequences in multi-fasta format\n";

    my $tempgp_fa = $work_dir . $species . $type . ".gp.fa";

    $tempgp_fa = TblToFasta( $tempgp_tbl, $tempgp_fa );

    unlink $pretblgp;

    print STDERR "\nSet up files for optimization\n\n";
    $my_command = "gawk '{print \$1,length(\$2)}' $tempgp_tbl | sort -k1,1 ";
    my $seqslenggp = capture($my_command);

    #  ` gawk '{print \$1,length(\$2)}' $tempgp_tbl | sort -k1,1 `;    ##XX

    my $tempseqlen = $work_dir . $species . $type . ".gp_cds_length";
    open( $fh_FOUT, ">", "$tempseqlen" ) or croak "Failed here";
    print $fh_FOUT "$seqslenggp";
    close $fh_FOUT;

    my $cdsgp = "";
    open( $fh_LOCID,
"./bin/gff2cds.awk source=\"annotations\" $tempgp_gff | sort -k1,1 | join $tempseqlen - |"
    );
    while (<$fh_LOCID>) {
        $cdsgp .= $_;
    }
    close $fh_LOCID;

    my $tempcdsgp = $work_dir . $species . $type . ".cds_gp";
    open( $fh_FOUT, ">", "$tempcdsgp" ) or croak "Failed here";
    print $fh_FOUT "$cdsgp";
    close $fh_FOUT;

    my $gffgpeval = "";
    open( $fh_LOCID,
"gawk 'BEGIN{while (getline<ARGV[1]>0){len[\$1]=\$2;};ARGV[1]=\"\";OFS=\"\\t\";}{if (NR==1) {ant=\$1;print \$1,\$2,\"Sequence\",1,len[ant],\"\.\",\"\.\",\"\.\",\"\.\"};if (\$1!=ant) {print \"\#\$\";ant=\$1;print \$1,\$2,\"Sequence\",1,len[ant],\"\.\",\"\.\",\"\.\",\"\.\"}; print }' $tempcdsgp $tempgp_gff |"
    );
    while (<$fh_LOCID>) {
        $gffgpeval .= $_;

    }
    close $fh_LOCID;

    my $tempevalgpgff = $work_dir . $species . $type . ".gp_eval_gff";
    open( $fh_FOUT, ">", "$tempevalgpgff" ) or croak "Failed here";
    print $fh_FOUT "$gffgpeval";
    close $fh_FOUT;

    if ($contigopt) {

        #my $tempgp_fa = $species.$type.".gp.fa";

        #$tempgp_fa = FastaToTbl($tempgp_tbl,$tempgp_fa);

        my @tabulargp = split( /\n/, $tblgp );
        my $seq = "";
        foreach my $line (@tabulargp) {
            chomp $line;
            my @f = split " ", $line;
            $seq .= $f[1];
        }
        my $lengp       = length($seq);
        my $foldedseqgp = fold4fasta($seq);
        my $tempfastagpcontig =
          $work_dir . $species . $type . ".combined.gp.fa";
        open( $fh_FOUT, ">", "$tempfastagpcontig" ) or croak "Failed here";
        print $fh_FOUT ">$species\n$foldedseqgp\n";
        close $fh_FOUT;

        my $temptabulargpcontig =
          $work_dir . $species . $type . ".combined.gp.tbl";
        open( $fh_FOUT, ">", "$temptabulargpcontig" ) or croak "Failed here";
        print $fh_FOUT "$species\t$seq\n";
        close $fh_FOUT;
        $my_command = "gawk '{print \$1,length(\$2)}' $temptabulargpcontig";
        my $seqslengcontiggp = capture($my_command);

        #` gawk '{print \$1,length(\$2)}' $temptabulargpcontig `;

        my $tempseqlencontig =
          $work_dir . $species . $type . ".gp_cds_contig_length";
        open( $fh_FOUT, ">", "$tempseqlencontig" ) or croak "Failed here";
        print $fh_FOUT "$seqslengcontiggp";
        close $fh_FOUT;

        my $gpcontig = "";
        open( $fh_LOCID,
"./bin/multiple_annot2one.awk species=$species leng=$lengp $tempcdsgp |"
        );
        while (<$fh_LOCID>) {

            $gpcontig .= $_;
        }
        close $fh_LOCID;

        my $tempgff2gpcontig = $work_dir . $species . $type . ".contig.gp.cds";

        open( $fh_FOUT, ">", "$tempgff2gpcontig" ) or croak "Failed here";
        print $fh_FOUT "$gpcontig";
        close $fh_FOUT;

        my $cds2gffcontig = "";
        open( $fh_LOCID,
"./bin/cds2gff.awk $tempgff2gpcontig | gawk 'BEGIN{OFS=\"\\t\";}{if (NR==1){print \"$species\",\"annotations\",\"Sequence\",\"1\",$lengp,\".\",\".\",\".\",\".\";print}else {print}}' - | "
        );
        while (<$fh_LOCID>) {
            $cds2gffcontig .= $_;
        }
        close $fh_LOCID;

        my $tempgp_cdsgff_contig_eval =
          $work_dir . $species . $type . ".cds_gp_contig.eval.gff";
        open( $fh_FOUT, ">", "$tempgp_cdsgff_contig_eval" )
          or croak "Failed here";
        print $fh_FOUT "$cds2gffcontig";
        close $fh_FOUT;

        return [
            $tempgp_cdsgff_contig_eval, $tempfastagpcontig,
            $temptabulargpcontig,       $tempseqlencontig
        ];

    }
    elsif ( !$contigopt ) {
        return [ $tempevalgpgff, $tempgp_fa, $tempgp_tbl, $tempseqlen ];
    }

}    #processSequences optimization

## GETGENES FUNCTION: EXTRACT FLANKED SEQUENCES FROM GENE MODELS FOR LATER OPTIMIZATION
sub GetGenes {

    my ( $path2gpath, $genesfname, $work_dir, $outtblgp ) = @_;

#print STDERR "IN FUNCTION: $path2gpath : $genesfname : $path : OUT: $outtblgp\n\n";

    my $nonred     = 0;
    my $onlynonred = 0;
    my $prevgene   = "x";
    my $prevchro   = "x";
    my $trail      = "";

    #my %genenames; unused var
    $outtblgp = "outputtbl.txt";

    chomp($path2gpath);
    chomp($genesfname);
    chomp($work_dir);
    chomp($outtblgp);

    #if ( !open( REFGENE, "< $genesfname" ) ) {
    #print "getgenes: impossible to open $genesfname\n";
    #exit(1);
    #}

    #if ( !open( OUT, "> $outtblgp" ) ) {
    #print "getgenes: impossible to create $outtblgp\n";
    #exit(1);
    #}

    open( my $fh_REFGENE,   "<", "$genesfname" ) or croak "Failed here";
    open( my $fh_OUT_tblgb, ">", "$outtblgp" )   or croak "Failed here";
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
        my $subseq    = "";
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
            $chrotmptbl = FastaToTbl( $path2gpath . $chro, $chrotmptbl );

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
            my $subseq = "";

            #my $sublen = 0;
            foreach my $line (@tabular) {
                chomp $line;
                my @f = split " ", $line;

                #print STDERR "$f[0]\n";
                $subseq .= $f[1];

                #$sublen += length($f[1]);
                #$countlines++;
                #if ($sublen >= ($numseqs * $k +1)){last;}
            }
######added

            if ( $le[1] < $txEn ) {

                my $newlen = $le[1];
                $genomic = substr( $subseq, $txSt, ( $newlen - $txSt ) );

            }
            elsif ( $le[1] >= $txEn ) {
                $genomic = substr( $subseq, $txSt, $txLe );

            }

            # my $genomic = `$call`;
            # my $genomicLe = length($genomic);
            my $genomicLe = length($genomic);
            my $cdseq     = "";

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

                            $utrs  = lc($utrs);
                            $cdseq = $cdseq
                              . "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n";
                        }

                        $iex   = $i + 1;
                        $cdseq = $cdseq
                          . "$name\t$chro\t$exTy\t$iex\t$exLe\t$seq\t$exon[$i]\t$exon[$i+$exoC]\n";

                        if ($utrA) {
                            my $iutr = $i + 1;
                            my $utrs =
                              substr( $genomic, $utrS + $cdsoffset, $utrL );

                            $utrs  = lc($utrs);
                            $cdseq = $cdseq
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
                            $utrs  = reverse($utrs);
                            $cdseq = "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n"
                              . $cdseq;
                        }

                        $iex = $exoC - $i;
                        $seq =~ tr/acgt/tgca/;
                        $seq = reverse($seq);
                        $cdseq =
"$name\t$chro\t$exTy\t$iex\t$exLe\t$seq\t$exon[$i+$exoC]\t$exon[$i]\n"
                          . $cdseq;

                        if ($utrA) {
                            my $iutr = $exoC - $i;
                            my $utrs =
                              substr( $genomic, $utrS + $cdsoffset, $utrL );

                            $utrs = lc($utrs);
                            $utrs =~ tr/acgt/tgca/;
                            $utrs  = reverse($utrs);
                            $cdseq = "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n"
                              . $cdseq;
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
                                $iIn   = $j + 1;
                                $cdseq = $cdseq
                                  . "$name\t$chro\tIntron\t$iIn\t$inLe\t$seq\n";
                            }
                            else {
                                $iIn = $exoC - $j - 1;
                                $seq =~ tr/acgt/tgca/;
                                $seq = reverse($seq);
                                $cdseq =
                                  "$name\t$chro\tIntron\t$iIn\t$inLe\t$seq\n"
                                  . $cdseq;
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
                        $cdseq =
                          $cdseq . "$name\t$chro\tUtr\t$iutr\t$exLe\t$utrs\n";
                    }
                    else {                   # reverse
                        my $iutr = $exoC - $i;

                        $utrs = lc($utrs);
                        $utrs =~ tr/acgt/tgca/;
                        $utrs = reverse($utrs);
                        $cdseq =
                          "$name\t$chro\tUtr\t$iutr\t$exLe\t$utrs\n" . $cdseq;
                    }

                }
            }

            print $fh_OUT_tblgb $cdseq;

        }
        elsif ($onlynonred) {
            print STDERR "$name\n";
        }

    }
    close $fh_OUT_tblgb;
    close $fh_REFGENE;

    return $outtblgp;

}    #getgenes XXX why does it not fold?

##GET BACKGROUND SEQUENCES (ex. 62 Kmer) used to compute log likelihoods

## GETKMATRIX FUNCTION (Splice site an Start codon PWMs)
sub getKmatrix {
    ## BUG do not relay on  jacknife/branch etc. use some: type ??
    #was my our?
    my (
        $true_seqs, $false_seqs, $order, $offset,
        $donor,     $accept,     $ATG,   $branch,
        $start,     $end,        $jacknife
    ) = @_;

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

    if ($branch) {
        $matrix_type = 'branch';
    }
    my $original_offset = $offset;
    my @prof            = ();
    my $tempinfolog;
    my @orders  = (qw(order-0 di tri order-4 order-5 order-6 order-7 order-8));
    my $ordname = $orders[$order];
    my $sort    = "sort -n";
    $sort = "sort -k1,1n -k2,2" if $order > 1;

    #	my @info = ($offset-1,$offset+1);
    my $prof_len      = 0;
    my $info_thresh   = "";           #bits
    my $pid           = $$;
    my $true_seq_name = $true_seqs;
    $true_seq_name =~ s/\.tbl$//;
    my $false_seq_name = $false_seqs;
    $false_seq_name =~ s/\.tbl$//;
## Open true sequences
    #	print STDERR "$true_seqs (true)\n";
    open( my $fh_true_seq, "<", "$true_seqs" ) or croak "Failed here";
    $_ = <$fh_true_seq>;
    my @t   = split;
    my $len = length( $t[1] );
    close $fh_true_seq;

## Open false (background???) sequences
    #	print STDERR "$false_seqs (false)\n";
    open( my $fh_FALSE_SEQ, "<", "$false_seqs" )
      or croak "Couldn't open $false_seqs: $!\n";
    $_ = <$fh_FALSE_SEQ>;
    my @f    = split;
    my $len2 = length( $f[1] );
    close $fh_FALSE_SEQ;

    #	die "$len != $len2\n" if $len != $len2;
    my $true_seq_freq_fn  = $work_dir . basename($true_seq_name) . ".freq";
    my $false_seq_freq_fn = $work_dir . basename($false_seq_name) . ".freq";

#my $subtracted_true_false_freq_fn = $work_dir . basename($true_seq_name) . "_" . basename($false_seq_name).freq_subtr";
    my $my_freq_subtract_fn =
        $work_dir
      . basename($true_seq_name) . "_"
      . basename($false_seq_name)
      . ".information";

    run("./bin/frequency.py 1 $true_seqs  >  $true_seq_freq_fn");
    run("./bin/frequency.py 1 $false_seqs >  $false_seq_freq_fn");

    #donor_gawk $1<=38 && \$1>=25
    #acceptor   $1<=33 && \$1>=2
    #start      $1<=37 && \$1>=25
    #branch     $1<=41 && \$1>=28'

    my $my_freq_field_limit_1 = 0;
    my $my_freq_field_limit_2 = 0;

    my %frequency_thresholds = (
        donor    => [ 38, 25 ],
        acceptor => [ 33, 2 ],
        ATG      => [ 37, 25 ],
        branch   => [ 41, 28 ],
    );

    $my_freq_field_limit_1 = $frequency_thresholds{$matrix_type}[0];
    $my_freq_field_limit_2 = $frequency_thresholds{$matrix_type}[1];

    #~ if ($donor) {
    #~ $my_freq_field_limit_1 = 38;
    #~ $my_freq_field_limit_2 = 25;
    #~ }
    #~ if ($accept) {
    #~ $my_freq_field_limit_1 = 33;
    #~ $my_freq_field_limit_2 = 2;    ### BUG !!! ???
    #~ }

    #~ if ($star) {
    #~ $my_freq_field_limit_1 = 37;
    #~ $my_freq_field_limit_2 = 25;
    #~ }

    #~ if ($branch) {
    #~ $my_freq_field_limit_1 = 41;
    #~ $my_freq_field_limit_2 = 28;
    #~ }
    my $my_command_A =
      "./bin/information.py  $true_seq_freq_fn $false_seq_freq_fn ";
    my $my_command_B =
"| gawk 'NF==2 && \$1<=$my_freq_field_limit_1 && \$1>=$my_freq_field_limit_2'"; 

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
    my $my_T_logratio_freq_matrix_fn =
      $work_dir . basename($true_seq_name) . ".log.$ordname.matrix";

    if ( !$order ) {
        $my_command =
"gawk -f ./bin/logratio_zero_order.awk $false_seq_freq_fn $true_seq_freq_fn > $my_T_logratio_freq_matrix_fn";
        say "\n $my_command \n";
        run($my_command);

    }
    else {
        $my_command =
"gawk -f ./bin/Getkmatrix.awk $order $len $true_seqs | $sort > $my_True_freq_matrix_fn";
        say "\n $my_command \n";
        run($my_command);

#~ run(
#~ " gawk -f ./bin/Getkmatrix.awk $order $len $true_seqs | $sort > $true_seq_name.$ordname-matrix"
#~ );

        $my_command =
"gawk -f ./bin/Getkmatrix.awk $order $len2 $false_seqs | $sort > $my_False_freq_matrix_fn ";
        say "\n $my_command \n";
        run($my_command);

        $my_command =
"gawk -f ./bin/logratio_kmatrix.awk $my_False_freq_matrix_fn $my_True_freq_matrix_fn > $my_T_logratio_freq_matrix_fn";
        say "\n $my_command \n";
        run($my_command);

    }

    #need to check output and then go on
## draw bit score bar graph function (nested, local)

    local *BitScoreGraph = sub {

        my ( $infooutput, $info_thresh, $offset ) = @_;
        my @info = ( $offset - 1, $offset + 1 );
        open( my $fh_INFO, "<", "$infooutput" ) or croak "Failed here";
        while (<$fh_INFO>) {
            next if m/^#/;
            last if m/^\s/;
            last if m/^[^\d]/;
            chomp;
            my @fields = split;
            printf STDERR "%2s %2.2f %s", ( $fields[0], $fields[1], "=" x int( $fields[1] * 30 ) );
            if ( $fields[1] > $info_thresh ) {
                push( @info, $fields[0] );
            }
            print STDERR "\n";
        }
        close $fh_INFO;
        print STDERR "\n BitScoreGraph \n";
        
        my @sortedinfo = sort numerically @info;
        my $start      = ( shift @sortedinfo );
        $start = 1 if $start < 1;
        my $end = pop @sortedinfo;

        return ( $start, $end );
    };    #end BitScoreGraph

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
    #~ if ( $donor && !$jacknife ) {

    #~ $info_thresh = "0.15";

    #~ ( $start, $end ) = BitScoreGraph( $tempinfolog, $info_thresh, $offset );

    #~ }
    #~ elsif ( $accept && !$jacknife ) {

    #~ $info_thresh = "0.04";

    #~ ( $start, $end ) = BitScoreGraph( $tempinfolog, $info_thresh, $offset );

    #~ }
    #~ elsif ( $star && !$jacknife ) {

    #~ $info_thresh = "0.15";

    #~ ( $start, $end ) = BitScoreGraph( $tempinfolog, $info_thresh, $offset );

    #~ }
    #~ elsif ( $branch && !$jacknife ) {

    #~ $info_thresh = "0.3";

    #~ ( $start, $end ) = BitScoreGraph( $tempinfolog, $info_thresh, $offset );

    #~ }

    if ( !$jacknife && $interactive ) {
        my $resp  = "";
        my $sline = "";
        my $eline = "";
        do {
            print STDERR
"Automatically selected subsequence from $start to $end to use in profile.\nDo you prefer to change the start or end? ";
            $resp = readline(STDIN);
        } while ( $resp !~ /^(yes|y)|(n|no)$/i );
        if ( $resp =~ /^(yes|y)/i ) {
            do {
                print STDERR "\nType new start: ";
                $sline = readline(STDIN);
            } while ( $sline !~ /(\d+)/ );
            $start = $1;
            do {
                print STDERR "\nType new end: ";
                $eline = readline(STDIN);
            } while ( $eline !~ /(\d+)/ || $eline <= $sline );
            $end = $1;
        }
    }    #if not jacknife

    $offset = $offset - $order;

    #$start = $start - $order;
    if ( !$jacknife ) {
        $end = $end - $order;
    }

    #	if ( $start < 1 ) {
    #	    $start = 1;
    #	}
    $offset = $offset - $start + 1;
    print STDERR "end:$end offset:$offset start:$start\n";
    if ( !$jacknife ) {
        print STDERR
          "new offset: $offset\nnew start: $start\nnew order: $order\n";
    }

    #my $my_Tlograt_summatrix_fn;
    #my $my_T_dimatrixdonor_4param_fn;
    my $my_Tlograt_summatrix_fn = $my_T_logratio_freq_matrix_fn . ".submatrix";
    my $my_T_dimatrixdonor_4param_fn =
      $my_T_logratio_freq_matrix_fn . ".dimatrix_4param";

    if ( $order >= 1 && $donor ) {

        my $preoffset = $offset + 2;
        my $newoffset = $offset + 3;
        my $posoffset = $offset + 4;

        #my_True_dimatrixdonor_4param_fn

        run(
" gawk -f ./bin/submatrix.awk $start $end $my_T_logratio_freq_matrix_fn > $my_Tlograt_summatrix_fn"
        );

        run(
"./bin/preparedimatrixdonor4parameter.awk $preoffset $newoffset $posoffset $my_Tlograt_summatrix_fn > $my_T_dimatrixdonor_4param_fn"
        );

# print STDERR "submatrix.awk $start $end $true_seq_name-log.$ordname-matrix/$true_seq_name-log-info.$ordname-matrix";

    }
    elsif ( $order >= 1 && $accept ) {

        my $preoffset = $offset - 1;
        my $newoffset = $offset;
        my $posoffset = $offset + 1;

        run(
" gawk -f ./bin/submatrix.awk $start $end $my_T_logratio_freq_matrix_fn > $my_Tlograt_summatrix_fn"
        );
        run(
" ./bin/preparedimatrixacceptor4parameter.awk $preoffset $newoffset $posoffset $my_Tlograt_summatrix_fn > $my_T_dimatrixdonor_4param_fn"
        );

#	  print STDERR "submatrix.awk $start $end $true_seq_name-log.$ordname-matrix/$my_True_dimatrixdonor_4param_fn";

    }
    elsif ( $order >= 2 && $ATG ) {

        my $preoffset = $offset - 2;
        my $newoffset = $offset - 1;
        my $posoffset = $offset;

        run(
" gawk -f ./bin/submatrix.awk $start $end $my_T_logratio_freq_matrix_fn > $my_Tlograt_summatrix_fn"
        );
        run(
" ./bin/preparetrimatrixstart4parameter.awk $preoffset $newoffset $posoffset $my_Tlograt_summatrix_fn > $my_T_dimatrixdonor_4param_fn"
        );
    }
    else {

# print STDERR "$path/submatrix_order0.awk $start $end $true_seq_name-log.$ordname-matrix\n";

        run(
" gawk -f ./bin/submatrix_order0.awk $start $end $my_T_logratio_freq_matrix_fn > $my_T_dimatrixdonor_4param_fn"
        );

    }
## CREATE DATA STRUCTURE CONTAINING MATRIX OF INTEREST

    open( my $fh_PROF, "<", "$my_T_dimatrixdonor_4param_fn" )
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

    for ( my $i = 0 ; $i < length($seq) ; $i += 60 ) {
        my $s = substr( $seq, $i, 60 );

        $foldedseq = $foldedseq . $s . "\n";

        #print STDERR $foldedseq;
    }

    return $foldedseq;
}

#Optimize parameter file
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

    open( $fh_SOUT, ">", "$species.OptimizeParameter.log" )
      or croak "Failed here";

    if ( !$branchswitch ) {
        print STDERR
          "\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\n\n";
        print $fh_SOUT
          "\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\n\n";
    }
    else {
        print STDERR
"\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\nMinBranch range : $iMin to $fMin\nAccCtx range : $iAccCtx to $fAccCtx\n\n";
        print $fh_SOUT
"\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\nMinBranch range : $iMin to $fMin\nAccCtx range : $iAccCtx to $fAccCtx\n\n";

    }

    for ( $IeWF = $IeWFini ; $IeWF <= $FeWF ; $IeWF += $deWF ) {
        print STDERR "eWF: $IeWF\noWF: ";

        for ( $IoWF = $IoWFini ; $IoWF <= $FoWF ; $IoWF += $doWF ) {
            if ( !$branchswitch ) {
                print STDERR "$IoWF  ";
            }
            if ($branchswitch) {
                print STDERR "$IoWF\nMinDist:  ";
                for ( $iMin = $iMinini ; $iMin <= $fMin ; $iMin += $dMin ) {
                    print STDERR "$iMin\nAccCtx:  ";

                    for (
                        $iAccCtx = $iAccCtxini ;
                        $iAccCtx <= $fAccCtx ;
                        $iAccCtx += $dAccCtx
                      )
                    {
                        print STDERR "$iAccCtx  ";
                        my $param = Geneid::Param->new();
                        $param->readParam("$species.geneid.param");

                        for ( my $i = 0 ; $i < $param->numIsocores ; $i++ ) {
                            if (
                                !defined @{ $param->isocores }[$i]
                                ->Exon_weights(
                                    [ $IeWF, $IeWF, $IeWF, $IeWF ]
                                )
                              )
                            {
                                croak "error in setting exon weights\n";
                            }
                            if ( !defined @{ $param->isocores }[$i]
                                ->Exon_factor( [ $IoWF, $IoWF, $IoWF, $IoWF ] )
                              )
                            {
# if (!defined @{$param->isocores}[$i]->Exon_factor([0.33,$bestIoWF,$bestIoWF,0.33])) {
                                croak "error in setting exon weights\n";
                            }
                            if (
                                !defined @{ $param->isocores }[$i]
                                ->Site_factor(
                                    [
                                        1 - $IoWF, 1 - $IoWF,
                                        1 - $IoWF, 1 - $IoWF
                                    ]
                                )
                              )
                            {
# if (!defined @{$param->isocores}[$i]->Site_factor([0.45,1-$bestIoWF,1-$bestIoWF,0.45])) {
                                croak "error in setting exon weights\n";
                            }
                            if (
                                !defined @{ $param->isocores }[$i]
                                ->set_profile(
                                    'Branch_point_profile', $prof_len_bra,
                                    $fxdbraoffset,          -50,
                                    0,                      0,
                                    1,                      $iAccCtx,
                                    $iMin,                  0,
                                    0,                      $branchmatrix
                                )
                              )
                            {
                                croak "error in setting profile\n";
                            }
                        }

                        $param->writeParam("$species.geneid.param.temp");
                        my $my_command;

                        #my $fh_geneid;
                        #my $fname_geneid;
                        print "\nbefore running my_command(s) geneid \n";
                        my $fh_geneid    = File::Temp->new();
                        my $fname_geneid = $fh_geneid->filename;
                        print "\ntemp geneid file: $fname_geneid \n";
                        $my_command =
"./bin/geneid -GP ${newparam}.temp $gpfa > $fname_geneid";
                        run($my_command);
                        my $fh_gawk_out    = File::Temp->new();
                        my $fname_gawk_out = $fh_gawk_out->filename;
                        $my_command =
"cat $fname_geneid | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' > $fname_gawk_out";
                        run($my_command);
                        $my_command =
"egrep -wv 'exon'  $fname_gawk_out > $tmp_dir/Predictions.${newparam}.gff";
                        run($my_command);

#` ./bin/geneid -GP ${newparam}.temp $gpfa | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' > $tmp_dir/Predictions.${newparam}.gff`;
                        print "\nafter running my_command(s) geneid \n";

                        my @evaluation_output = split " ",
` ./bin/evaluation -sta $tmp_dir/Predictions.${newparam}.gff $gpgff | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $IoWF, $IeWF, $iAccCtx, $iMin, \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}'`;

                        push( @evaluation_total, \@evaluation_output );

                    }    #$iAccCtx for loop
                    unless ( $iMin == $fMin ) {
                        print STDERR "\nMinDist:  ";
                    }
                    $iAccCtx = $iAccCtxini;

                }    #$iMin for loop
                unless ( $IoWF == $FoWF ) {
                    print STDERR "\noWF:  ";
                }
                $iAccCtx = $iAccCtxini;
                $iMin    = $iMinini;

            }    #if $branchswitch
            elsif ( !$branchswitch ) {

                my $param = Geneid::Param->new();
                $param->readParam("$species.geneid.param");

                for ( my $i = 0 ; $i < $param->numIsocores ; $i++ ) {
                    if ( !defined @{ $param->isocores }[$i]
                        ->Exon_weights( [ $IeWF, $IeWF, $IeWF, $IeWF ] ) )
                    {
                        croak "error in setting exon weights\n";
                    }
                    if ( !defined @{ $param->isocores }[$i]
                        ->Exon_factor( [ $IoWF, $IoWF, $IoWF, $IoWF ] ) )
                    {
 #   if (!defined @{$param->isocores}[$i]->Exon_factor([0.4,$IoWF,$IoWF,0.4])) {
                        croak "error in setting exon weights\n";
                    }
                    if (
                        !defined @{ $param->isocores }[$i]->Site_factor(
                            [ 1 - $IoWF, 1 - $IoWF, 1 - $IoWF, 1 - $IoWF ]
                        )
                      )
                    {
#	  if (!defined @{$param->isocores}[$i]->Site_factor([0.55,1-$IoWF,1-$IoWF,0.55])) {
                        croak "error in setting exon weights\n";
                    }

                }

                $param->writeParam("$species.geneid.param.temp");
                ###
                my $fh_geneid    = File::Temp->new();
                my $fname_geneid = $fh_geneid->filename;
                print "\ntemp geneid file: $fname_geneid \n";
                my $my_command =
                  "./bin/geneid -GP ${newparam}.temp $gpfa > $fname_geneid";
                run($my_command);
                $my_command =
"cat $fname_geneid | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' > $tmp_dir/Predictions.${newparam}.gff";
                run($my_command);

#` ./bin/geneid -GP ${newparam}.temp $gpfa | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' > $tmp_dir/Predictions.${newparam}.gff`;
## BUG very complex comand line
                my @evaluation_output = split " ",
` ./bin/evaluation -sta $tmp_dir/Predictions.${newparam}.gff $gpgff | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $IoWF, $IeWF, \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}' `;

                push( @evaluation_total, \@evaluation_output );

            }    #elseif

        }

        # $iAccCtx = $iAccCtxini;
        # $iMin=$iMinini;
        $IoWF = $IoWFini;
        print STDERR "\n";

    }

    return \@evaluation_total;

}    # end sub optimize parameter file

sub BuildOptimizedParameterFile {

    my ( $evalarray, $branchswitch, $prof_len_bra,
        $fxdbraoffset, $branchmatrix ) = @_;
    my @sortedeval     = ();
    my @evaluationinit = ();
    my $bestIoWF       = "";
    my $bestIeWF       = "";
    my $bestMin        = "";
    my $bestAcc        = "";
    my $fh_SOUT;
    open( $fh_SOUT, ">", "$species.BuildOptimizedParameterFile.log" )
      or croak "Failed here";
    print $fh_SOUT "\n input to func
             \n1:"
      . $evalarray . "\n2 "
      . $branchswitch . "\n3 "
      . $prof_len_bra . "\n4 "
      . $fxdbraoffset . "\n5 "
      . $branchmatrix . "\n";

    if ( !$branchswitch ) {
        ## BUG ???
        @sortedeval = sort sorteval @$evalarray;

        $bestIoWF = $sortedeval[0][0];    #0.2
        $bestIeWF = $sortedeval[0][1];    #-3.5

        print STDERR "\nBest performance obtained using IoWF: "
          . $sortedeval[0][0]
          . " and IeWF: "
          . $sortedeval[0][1] . "\n";
        print $fh_SOUT "\nBest parameter file performance obtained using oWF: "
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

            print $fh_SOUT join( "\t", @$eval_ref ), "\n";

        }

## FOUR BEST PERFORMANCE RESULTS AFTER OPTIMIZATION
        for ( my $i = 0 ; $i <= 2 ; $i++ ) {
            print STDERR join( "\t", @{ $sortedeval[$i] } ), "\n";
        }
############

## BUILD BEST PERFORMANCE (OPTIMIZED) PARAMETER FILE (USE VALUES OBTAINED ABOVE)

        my $param = Geneid::Param->new();
        $param->readParam("$species.geneid.param");

        for ( my $i = 0 ; $i < $param->numIsocores ; $i++ ) {
            if (
                !defined @{ $param->isocores }[$i]->Exon_weights(
                    [ $bestIeWF, $bestIeWF, $bestIeWF, $bestIeWF ]
                )
              )
            {
                croak "error in setting exon weights\n";
            }
            if ( !defined @{ $param->isocores }[$i]
                ->Exon_factor( [ $bestIoWF, $bestIoWF, $bestIoWF, $bestIoWF ] )
              )
            {
#	 if (!defined @{$param->isocores}[$i]->Exon_factor([0.4,$bestIoWF,$bestIoWF,0.4])) {
                croak "error in setting exon weights\n";
            }
            if (
                !defined @{ $param->isocores }[$i]->Site_factor(
                    [
                        1 - $bestIoWF,
                        1 - $bestIoWF,
                        1 - $bestIoWF,
                        1 - $bestIoWF
                    ]
                )
              )
            {
#	 if (!defined @{$param->isocores}[$i]->Site_factor([0.55,1-$bestIoWF,1-$bestIoWF,0.55])) {
                croak "error in setting exon weights\n";
            }
        }

        #write new parameter file (optimized)
        $param->writeParam("$species.geneid.optimized.param");

        print STDERR
"\n\nNew optimized parameter file named: $species.geneid.optimized.param \n";
        print $fh_SOUT
"\n\nNew optimized parameter file named: $species.geneid.optimized.param \n";

        return [ $bestIeWF, $bestIoWF, 0, 0, \@evaluationinit ];

    }
    elsif ($branchswitch) {

        @sortedeval = sort sortevalbranch @$evalarray;

        $bestIoWF = $sortedeval[0][0];    #0.2
        $bestIeWF = $sortedeval[0][1];    #-3.5
        $bestAcc  = $sortedeval[0][2];
        $bestMin  = $sortedeval[0][3];

        print STDERR "\nBest performance obtained using IoWF: "
          . $sortedeval[0][0]
          . " , IeWF: "
          . $sortedeval[0][1]
          . ", MinimalBranchDist: "
          . $bestMin
          . ", Optimal Branch Context: "
          . $bestAcc . "\n";
        print $fh_SOUT "best performance obtained using IoWF: "
          . $sortedeval[0][0]
          . " , IeWF: "
          . $sortedeval[0][1]
          . ", MinimalBranchDist: "
          . $bestMin
          . ", Optimal Branch Context: "
          . $bestAcc . "\n";

        #INITIALIZE ARRAY WITH EVALUATION PARAMETERS
        @evaluationinit =
          (
            qw(oWF eWF AccCtx MinD SN SP CC SNe SPe SNSP SNg SPg SNSPg raME raWE)
          );

        print STDERR
"(Sorted performance results (Three best performance estimates) for different values of oWF, eWF, AccCtx and MinD)\n\n"
          . join( "\t", @evaluationinit ), "\n";
        print $fh_SOUT
"(Sorted performance results (best to worst) for different values of oWF, eWF, AccCtx and MinD)\n\n"
          . join( "\t", @evaluationinit ), "\n";

        foreach my $eval_ref (@sortedeval) {

            print $fh_SOUT join( "\t", @$eval_ref ), "\n";

        }

###THREE BEST PERFORMANCE RESULTS AFTER OPTIMIZATION
        for ( my $i = 0 ; $i <= 2 ; $i++ ) {
            print STDERR join( "\t", @{ $sortedeval[$i] } ), "\n";
        }
############

        ## BUILD BEST PERFORMANCE (OPTIMIZED) PARAMETER FILE (USE VALUES OBTAINED ABOVE)

        my $param = Geneid::Param->new();
        $param->readParam("$species.geneid.param");

        for ( my $i = 0 ; $i < $param->numIsocores ; $i++ ) {
            if (
                !defined @{ $param->isocores }[$i]->Exon_weights(
                    [ $bestIeWF, $bestIeWF, $bestIeWF, $bestIeWF ]
                )
              )
            {
                croak "error in setting exon weights\n";
            }
            if ( !defined @{ $param->isocores }[$i]
                ->Exon_factor( [ $bestIoWF, $bestIoWF, $bestIoWF, $bestIoWF ] )
              )
            {
#   if (!defined @{$param->isocores}[$i]->Exon_factor([0.4,$bestIoWF,$bestIoWF,0.4])) {
                croak "error in setting exon weights\n";
            }
            if (
                !defined @{ $param->isocores }[$i]->Site_factor(
                    [
                        1 - $bestIoWF,
                        1 - $bestIoWF,
                        1 - $bestIoWF,
                        1 - $bestIoWF
                    ]
                )
              )
            {
#  if (!defined @{$param->isocores}[$i]->Site_factor([0.55,1-$bestIoWF,1-$bestIoWF,0.55])) {
                croak "error in setting exon weights\n";
            }
            if (
                !defined @{ $param->isocores }[$i]->set_profile(
                    'Branch_point_profile', $prof_len_bra, $fxdbraoffset,
                    -50, 0, 0, 1, $bestAcc, $bestMin, 0, 0, $branchmatrix
                )
              )
            {
                croak "error in setting profile\n";
            }
        }

        #write new parameter file (optimized)
        $param->writeParam("$species.geneid.optimized.param");

        print STDERR
"\nNew optimized parameter file named: $species.geneid.optimized.param \n";
        print $fh_SOUT
"\nNew optimized parameter file named: $species.geneid.optimized.param \n";

        return [ $bestIeWF, $bestIoWF, $bestAcc, $bestMin, \@evaluationinit ];

    }    #if branch switch
    close $fh_SOUT;
    return 1;
}

sub EvaluateParameter {

    my ( $gpfa, $gpgff, $newparam ) = @_;

` ./bin/geneid -GP $newparam $gpfa | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' | egrep -wv 'exon' > test.predictions.${newparam}.gff `;

    my @evaluation_test = split " ",
` ./bin/evaluation -sta test.predictions.${newparam}.gff $gpgff | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}' `;

    return \@evaluation_test;

}    # evaluate parameter function

sub WriteStatsFile {

    my (
        $species,        $sout,           $outintron,      $outcds,
        $outgff,         $inframe,        $inframeeval,    $seqsused,
        $totalnoncandon, $totalnoncanacc, $totalnoncansta, $markovmodel,
        $totalcoding,    $totalnoncoding, $stdo,           $endo,
        $stac,           $enac,           $stst,           $enst,
        $stbr,           $enbr,           $branchsw,       $useallseqs
    ) = @_;
    my $my_command = "";
    
## OBTAIN GENE MODEL SET STATISTICS
## Open gene model object

    $param->geneModel( Geneid::GeneModel->new() );
    $param->geneModel->useDefault;
###########
    my $fh_SOUT;
    open( $fh_SOUT, ">", "test.WriteStatsFileLog.txt" ) or croak "Failed here";

    my $avgintron = "";
    my $sdintron  = "";

    my @intronlength = ` gawk '{print length(\$2)}' $outintron `;
    my ( $mean, $st ) = average( \@intronlength );

    #print STDERR "INTRON: $mean, $st\n";
    $my_command = "gawk '{print length(\$2)}' $outintron | sort -n";
    my @intronlist = capture($my_command);
    #my @intronlist = ` gawk '{print length(\$2)}' $outintron | sort -n `;

    my $totintrons = scalar(@intronlist);
    my @intronlen  = ();
    my $intr       = "";
    for ( my $i = 0 ; $i <= scalar(@intronlist) - 1 ; $i++ ) {

        $intr = $intronlist[$i];
        chomp $intr;
        push( @intronlen, $intr );
    }

    my @slice1 = @intronlen[ 0 .. 5 ];
    my @slice2 =
      @intronlen[ ( scalar(@intronlen) - 5 ) .. ( scalar(@intronlen) - 1 ) ];

    my $shortintron = $intronlist[0] - $intronlist[0] * (0.25);
    chomp $shortintron;
    if ( $shortintron > 40 ) { $shortintron = 40; }

    #my $longintron =  $intronlist[$totintrons - 1];
    my $longintron =
      $mean + ( $st * 3 ) > 100000 ? 100000 : $mean + ( $st * 3 );
    chomp $longintron;
    if ($interactive) {
        my $answer = "";
        do {
            print
"\nDo you want to use this automatically selected range of introns ($shortintron to $longintron) in the geneid gene model (yes/no)?\n(Note that the 5 smallest introns were found to be: "
              . join( ", ", @slice1 )
              . " nucleotides long and the 5 longest introns: "
              . join( ", ", @slice2 )
              . " bases in length)\n";
            $answer = <ARGV>;
            chomp $answer;
        } while ( $answer !~ /^(yes|y)|(n|no)$/i );

        if ( $answer =~ /^(no|n)$/i ) {
            my $sline = "";
            my $eline = "";
            do {
                print STDERR "\nType new intron minimum length boundary: ";
                $sline = readline(STDIN);
            } while ( $sline !~ /(\d+)/ );
            $shortintron = $1;
            do {
                print STDERR "\nType maximum intron length boundary: ";
                $eline = readline(STDIN);
            } while ( $eline !~ /(\d+)/ || $eline <= $sline );
            $longintron = $1;
        }

    }
    my $minintergenic = 200;
    my $maxintergenic = 'Infinity';
    if ($interactive) {
        do {
            print
"\nDo you want to use this automatically selected intergenic distance range ($minintergenic to $maxintergenic) in the geneid gene model (yes/no)?\n";
            $answer = <ARGV>;
            chomp $answer;
        } while ( $answer !~ /^(yes|y)|(n|no)$/i );

        if ( $answer =~ /^(no|n)$/i ) {
            my $sline = "";
            my $eline = "";
            do {
                print STDERR "\nType new minimum intergenic distance: ";
                $sline = readline(STDIN);
            } while ( $sline !~ /(\d+)/ );
            $minintergenic = $1;
            do {
                print STDERR
"\nType maximum intergenic distance (type a number or 'Infinity'): ";
                $eline = readline(STDIN);
            } while ( $eline !~ /(\d+|Infinity)/ || $eline <= $sline );
            $maxintergenic = $1;

        }
    }
## use shortest and longest intron lengths in gene model of parameter file
    $param->geneModel->intronRange( $shortintron, $longintron );
    $param->geneModel->intergenicRange( $minintergenic, $maxintergenic );
###############################
    
    $my_command = "gawk '{print gsub(/[GC]/,\".\",\$2)/length(\$2)} $outcds";
    my @CDSGCcontent = capture($my_command);
        
#    my @CDSGCcontent =
#    `gawk '{print gsub(/[GC]/,".",\$2)/length(\$2)}' $outcds `;

    #print STDERR "@CDSGCcontent\n";
    my ( $meangc, $stgc ) = average( \@CDSGCcontent );

    #print STDERR "CDS: $meangc $stgc $outintron\n";
    $my_command = "gawk '{print gsub(/[GC]/,\".\",\$2)/length(\$2)} $outintron ";
    my @intronGCcontent = capture($my_command);
    
    #my @intronGCcontent =
    #  ` gawk '{print gsub(/[GC]/,".",\$2)/length(\$2)}' $outintron `;

    #print STDERR "@intronGCcontent\n";
    my ( $meangci, $stgci ) = average( \@intronGCcontent );

    #print STDERR "intron: $meangci $stgci\n";
    #BUG?
    #my $totexons = ` gawk '{print \$9}' $outgff | wc -l | gawk '{print \$1}' `;
    $my_command = "gawk '{print \$9}' $outgff | wc -l ";
    my $totexons = ` gawk '{print \$9}' $outgff | wc -l `;

    chomp $totexons;
    $totexons = int($totexons);
    my @exonspergene;
    @exonspergene =
      ` gawk '{print \$9}' $outgff | sort | uniq -c | gawk '{print \$1}' `;

    my ( $avgex, $stex ) = average( \@exonspergene );

    my @exonlength =
      ` egrep -v 'Single' $outgff | gawk '{len=\$5-\$4;print len}' - | sort `;

    my ( $avgle, $stle ) = average( \@exonlength );

    my $singlegenes = `egrep -c '(Single)' $outgff `;
    chomp $singlegenes;

    #print $fh_SOUT "GENE MODEL STATISTICS FOR $species\n\n";

    print $fh_SOUT
"\nA subset of $totalseqs4training sequences (randomly chosen from the $total_seqs gene models) was used for training\n\n";

    if ( !$useallseqs ) {
        print $fh_SOUT
"The user has selected to use $seqsused gene models (80 % of total) for training and to set aside $gffseqseval annotations (20 % of total) for evaluation\n\n";
    }
    else {
        print $fh_SOUT
"$total_seqs gene models were used for both training and evaluation\n\n";
    }

    if ( !$useallseqs ) {
        print $fh_SOUT
"$inframe of the gene models translate into proteins with in-frame stops within the training set and $inframeeval in the evaluation set (seqs removed).\n\n";
    }
    else {
        print $fh_SOUT
"$inframe of the gene models translate into proteins with in-frame stops within the training set.\n\n";
    }
    print $fh_SOUT
"There are $totalnoncandon non-canonical donors as part of the training set\n\n";
    print $fh_SOUT
"There are $totalnoncanacc non-canonical acceptors as part of the training set\n\n";
    print $fh_SOUT
"There are $totalnoncansta non-canonical start sites as part of the training set\n\n";
    print $fh_SOUT
"These gene models correspond to $totalcoding coding bases and $totalnoncoding non-coding bases\n\n";
    print $fh_SOUT
"Deriving a markov model for the coding potential of order $markovmodel\n\n";
    print $fh_SOUT
"The intronic sequences extracted from the gene models have an average length of $mean, with $st of SD\n";
    print $fh_SOUT
"Geneid can predict gene models having introns with a minimum length of $shortintron nucleotides and a maximum of $longintron bases (boundaries used in gene model) \n\n";
    print $fh_SOUT
"The minimum (user selected) intergenic distance was set to $minintergenic nucleotides whereas the maximum was set to $maxintergenic (boundaries used in gene model) \n\n";
    print $fh_SOUT
"The GC content of the exonic and intronic sequences is $meangc (SD $stgc) and $meangci (SD $stgci) respectively \n\n";
    print $fh_SOUT
      "The gene models used for training contain $totexons exons \n\n";
    print $fh_SOUT
      "The gene models average $avgex exons per gene (SD $stex)\n\n";
    print $fh_SOUT
"The average length of the exons (non-single) in the training set gene models is $avgle (SD $stle)\n\n";

    if ( !$useallseqs ) {
        print $fh_SOUT
"The training set includes $singlegenes single-exon genes (out of $seqsused ) gene models\n\n";
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
    return ( $shortintron, $longintron, $minintergenic, $maxintergenic );

}

sub average {
    my ($sequences) = @_;
    my $sum         = 0;
    my $total       = 0;
    my ( $mean, $st );

    foreach my $seq (@$sequences) {
        $sum += $seq;
        $total++;
    }

    $mean = $sum / $total;
    $mean = sprintf( "%.3f", $mean );

    $sum = 0;

    foreach my $seq (@$sequences) {
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
        chomp $_;
        my ( $n, $s ) = split( /\s+/, $_ );
        my ( $i, $e ) = ( 1, length($s) );
        print $fh_FOUT ">$n\n";
        while ( $i <= $e ) {
            print $fh_FOUT substr( $s, $i - 1, 60 ) . "\n";
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
        print $fh_FOUT_fasta ">$seq_name\n";
        print STDERR "#";
        while ( $base_number <= $seq_length ) {
            print $fh_FOUT_fasta substr( $seq, $base_number - 1, 60 ) . "\n";
            $base_number += 60;
        }
        close $fh_FOUT_fasta;
    }

    close $fh_IN_tbl;

    #close $fh_FOUT_fasta;
    return 1;
}

sub FastaToTbl {
    my ( $fa, $tblout ) = @_;

    open( my $fh_IN,   "<", "$fa" );
    open( my $fh_TOUT, ">", "$tblout" );

    #print STDERR "$fa INSIDE LOOP\n\n";
    my $count = 0;
    while (<$fh_IN>) {
        chomp;
        $_ =~ s/\|//;
        if ( $_ =~ /\>(\S+)/ ) {
            print $fh_TOUT "\n" if $count > 0;
            print $fh_TOUT $1 . "\t";
            $count++;
        }
        else {
            print $fh_TOUT $_;
        }
    }
    print $fh_TOUT "\n";

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

sub Translate {

    my ( $geneticcode_fn, $cds_fn, $outprot ) = @_;
    print STDERR "$geneticcode_fn in loop\n";
    my $frame = 0;

    #~ if ( !open(my $fh_FILEIN, "<", "$geneticcode" ) ) {
    #~ print "Translate: impossible to open genetic.code\n";
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
    #print "Translate: impossible to open $cds\n";
    #exit(1);
    #}

    open( my $fh_CDS_IN, "<", "$cds_fn" )  or croak "err:  $cds_fn";
    open( my $fh_POUT,   ">", "$outprot" ) or croak "Failed here";
    print STDERR "translating: $cds_fn \n";

    while (<$fh_CDS_IN>) {

        my $line = $_;

        my ( $name, $seq ) = split( /\s+/, $line );

        my $lseq = length($seq);

        print $fh_POUT "$name ";

        for ( my $i = $frame ; $i < $lseq - 1 ; $i += 3 ) {
            my $triplet = substr( $seq, $i, 3 );
            if ( $gencodeh{$triplet} ) {
                print $fh_POUT "$gencodeh{$triplet}";
            }
            else {
                print $fh_POUT "X";
            }

        }    #for
        print $fh_POUT "\n";
    }

    close $fh_CDS_IN;
    close $fh_POUT;
    say "\noutprot : $outprot \n";
    return $outprot;

}

## CONVERT GFF2 TO GENEID GFF
sub generalGFFtoGFFgeneid {

    my ( $gff, $species, $type ) = @_;
    my %G;
    my @G = ();

    my $geneidgff = $work_dir . $species . ${type} . ".geneid_gff";

    open( my $fh_GFF,    "<", "$gff" )       or croak "Failed here";
    open( my $fh_GFFOUT, ">", "$geneidgff" ) or croak "Failed here";
    while (<$fh_GFF>) {
        my ( $c, @f, $id );
        $c = ":";
        $_ =~ s/\|//;
        chomp;
        @f = split /\s+/o, $_;
        $id = $f[8];    #seq name i.e. 7000000188934730
        ( exists( $G{$id} ) ) || do {
            $c = "#";
            $G{$id} = [ @f[ 8, 0, 6 ], 0, @f[ 3, 4 ], [] ]
              ;    # [7000000188934730 1.4 - 0 46549	46680 ]
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
        ( $id, $ctg, $str, $nex, $go, $ge, $coords ) = @$GN;

       # print STDERR Data::Dumper->Dump([ \@$coords ], [ qw/ *coords / ])."\n";
        @coords = map { $_->[0], $_->[1] }
          sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
          map { $_ } @$coords;

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
            print $fh_GFFOUT join( "\t",
                $ctg, "$species", $feat, $coords[$ce], $coords[ ( $ce + 1 ) ],
                ".", $str, ".", $id )
              . "\n";
            $ce += 2;
            $c++;
        }
    }

    my $geneidgffsorted = "";
    open( my $fh_LOCID, "sort -s -k8,9 -k4,5n $geneidgff |" );
    while (<$fh_LOCID>) {
        $geneidgffsorted .= $_;
    }
    close $fh_LOCID;

####
    my $tempgeneidgffsorted =
      $work_dir . $species . $type . ".geneid.gff_sorted";
    open( my $fh_FOUT, ">", "$tempgeneidgffsorted" ) or croak "Failed here";
    print $fh_FOUT "$geneidgffsorted";
    close $fh_FOUT;

    close $fh_GFF;
    close $fh_GFFOUT;

    return $tempgeneidgffsorted;

    #exit(0);

}

## NO need to clean if we need to look at temp files to debug
#~ sub go_to_die() {
#~ ( warn "@_\n" ); #&& &clean_tmp();
#~ &clean_ext();    #unless $ext_flg;
#~ exit(1);
#~ }

## EVAL OPTIMIZATION SORTING FUNCTION

## NO need to clean if we need to look at temp files to debug
#~ sub clean_ext() {

#~ # Obtaining the list of extended files
#~ my @files = (qw( "" ));

#~ # Unlinking the temporary files if they exist
#~ foreach my $file (@files) {
#~ unlink $file if ( -e $file );
#~ }

#~ # rmdir "geneid_params" if (-e "geneid_params");

#~ }

## no need to look for patterns if all temp files are outside of the /tmp
#~ sub clean_tmp() {

#~ # Obtaining the list of temporary files
#~ opendir( my $dh_DIR, "$tmp_dir" ) or croak "err: $tmp_dir\n";
#~ my @files = map { "$tmp_dir/$_" } grep { /^$TMPROOT/ } readdir($dh_DIR);
#~ closedir($dh_DIR);

#~ foreach my $file (@files) {
#~ unlink $file if ( -e $file );
#~ }
#~ }

sub sorteval {

         $b->[7] <=> $a->[7]
      || $b->[10] <=> $a->[10]
      || $b->[4] <=> $a->[4]
      || $a->[11] <=> $b->[11]
      || $a->[12] <=> $b->[12]

}

sub sortevalbranch {

         $b->[9] <=> $a->[9]
      || $b->[12] <=> $a->[12]
      || $b->[6] <=> $a->[6]
      || $a->[13] <=> $b->[13]
      || $a->[14] <=> $b->[14]

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
    }

    #BASH AWK
    #system("which gff2ps > /dev/null;")
    #    && &go_to_die("The gff2ps package is not found or is not executable");
    say "/nNeccessary binaries are executable\n";
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
      or croak "Could not open file $input_tbl_fn' $!";

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
          or croak "Could not open file '$tmp_out_fn' $!";
        print $fh_out "$name $lengfasta\n";
        close $fh_out;

    }
    close $fh_input;
    return 1;
}    #end write_sizes_from_tbl_fn

## extracted from main flow
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
}    ###IF THERE IS A MEME-DISCOVERED BRANCH POINT PROFILE END

sub start_reduced {
#############################################################
## reduced/short version training starting with PWMs PLUS BACKGROUND#
#############################################################
    croak "no reduced training option in the simplified script";
    return 1;
}    #end reduced

sub getBackground {

    my ( $k, $fasta, $tbl, $numseqs, $background ) = @_;

    #my $fasta = $f;
    #my $tbl = "";
    my $countlines = 0;
    my $totalseqs  = num_of_lines_in_file($fasta);

    #my $totalseqs  = `egrep -c \"^>\" $fasta`;
    #chomp $totalseqs;

    open( my $fh_IN, "<", "$tbl" ) or croak "Failed here";
    my @tabular = ();
    while (<$fh_IN>) {
        push @tabular, $_;
    }
    close $fh_IN;

    print STDERR
      "\nThe total number of genomic sequences (in $fasta) is $totalseqs\n";

    if ( $totalseqs >= 1 ) {

        #print STDERR "in totalseqs if";
        my $seq    = "";
        my $sublen = 0;
        foreach my $line (@tabular) {
            chomp $line;
            my @f = split " ", $line;
            $seq .= $f[1];
            $sublen += length( $f[1] );
            $countlines++;
            if ( $sublen >= ( $numseqs * $k + 1 ) ) {
                last;
            }
        }

        #print STDERR "$seq";
        my $len = length($seq);

        # chomp $totalseqs;
        print STDERR
"The length of the new fasta file is $len\n(concatenation of $countlines fastas (out of $totalseqs))\n";
        open( my $fh_BACKGRND, ">", "$background" ) or croak "Failed here";
        my $row = 1;
        print STDERR
          "\nObtain background ($numseqs seqs) from fasta of $k nucleotides \n";
        for ( my $n = 1 ; $n <= $numseqs ; $n++ ) {
            my $kmer = "N";
            while ( $kmer =~ m/[^ACGTacgt]/ ) {
                my $rand = int( rand( $len - $k ) + 1 );
                $kmer = substr( $seq, $rand, $k );
            }

            #print STDERR  $n."..";
            print $fh_BACKGRND $row++ . "\t$kmer\n";
        }
        print STDERR "FINISHED OBTAINING BACKGROUND INFO\n";
        close $fh_BACKGRND;

    }

    #return $background;

}    # END getbackground function

