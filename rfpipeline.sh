#!/bin/bash

#python protparams.py $1
echo "Extracting features from $1 ..."
IFS=.
set $1
mkdir $1
echo "All feature data will be stored in folder '$1'"
#tr -d '{}[]()'\''":' < protparams.csv > protparamsdelim.csv
#Rscript frequency.r
#cp frequency.csv protparamsdelim.csv $1
#rm *.csv
cp ${1}.fasta $1
cd $1
mkdir ref
csplit -z ${1}.fasta '/^>/' '{*}' --suffix="%02d.fasta" --prefix=$1- -s
mv ${1}.fasta ref
cd ..
cp psipred/multiseqpsipred.sh $1
cp psipred/runpsipredplus $1
cp psipred/psipred.pl $1
cp psipred/psiparse.r $1
cd $1
pwd 
./multiseqpsipred.sh
Rscript psiparse.r
#Rscript dataset.r

