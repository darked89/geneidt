## install required Perl modules

i.e. by running:

`cpanm Carp::Always`

`cpanm  Carp::Assert::More`

`cpanm Data::Dumper`

`cpanm Cwd`

`cpanm Getopt::Long`

`cpanm File::Path`

`cpanm File::Temp`

`cpanm IPC::System::Simple`

##  compile C programs (make)

`src/c/evaluation
src/c/ssgff
src/c/pictogram
`
compile geneid (not included in this distribution)

## set the PERL5LIBS:
i.e. in fish shell
`set -gx PERL5LIB $HOME/perl5/lib/perl5/ $PERL5LIB`

## run pmarinus totalseqs4training
``./run_test_train.sh`

## clean up all training results
`clean_test_train_results.sh`
