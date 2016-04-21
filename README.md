# README #

To install the geneid_training from Bitbucket

### What is the status ###
Latest commit works  (80% training, 20% for eveluation mode)



* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* install required Perl modules
easiest is to use cpanm to get modules located in:
`$HOME/perl5/lib/perl5/`
 

`cpanm Carp`

`cpanm Carp::Always`

`cpanm Carp::Assert` 

`cpanm Data::Dumper`

`cpanm Cwd`

`cpanm Getopt::Long`

`cpanm File::Path`

`cpanm File::Basename`

`cpanm File::Temp` 

`cpanm IPC::System::Simple` 

`cpanm Readonly`

set the PERL5LIBS:
i.e. in fish shell
`set -gx PERL5LIB $HOME/perl5/lib/perl5/ $PERL5LIB`

## run pmarinus totalseqs4training
``./run_test_train.sh`

## clean up all training results
`clean_test_train_results.sh`

* C programs 
```
cd src/c/evaluation
make
cp -i ./bin/* ../../../bin/
cd ../../../

cd src/c/ssgff
make
cp -i ./bin/* ../../../bin/
cd ../../../

cd src/c/pictogram
make
cp -i pictogram ../../../bin/
cd ../../../
```
compile geneid (not included in this distribution)

### How do I run the test ###
`./run_test_train.sh`

### clean up all training results ###
`clean_test_train_results.sh`